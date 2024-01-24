#include <fstream>
#include <iostream>
#include <string>
#include <array>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

// Boilerplate common across inputs for distributed RIMP2 tests
json template_output() {
  json output = json::object();
  output["model"] = {
    {"method", "rimp2"},
    {"spin_configuration", "restricted"},
    {"fragmentation", true},
    {"basis", "cc-pVDZ"},
    {"aux_basis", "cc-pVDZ-RIFIT"}
  };

  output["system"] = {
    {"max_gpu_memory_mb", 30000}
  };

  output["keywords"] = json::object();
  output["keywords"]["scf"] = {
    {"niter", 50},
    {"ndiis", 8},
    {"scf_conv", 1e-08},
    {"convergence_metric", "energy"}
  };

  output["keywords"]["frag"] = {
    {"fragmentation_level", "trimer"},
    {"fragmented_energy_type", "total_energy"},
    {"ngpus_per_node", 4},
    {"cutoffs", {
      {"dimer", 40},
      {"trimer", 30}}
    }
  };

  output["keywords"]["guess"] = {
    {"superposition_monomer_densities", false}
  };

  output["driver"] = "energy";
  return output;
}

struct Atom {
  std::string symbol;
  std::array<double,3> coord;
};

struct Fragment {
  int charge;
  std::vector<Atom> atoms;
};

std::vector<Atom> parse_atoms(json molecule) {
  std::vector<std::string> symbols = molecule["symbols"];
  std::vector<double> coords = molecule["geometry"];
  std::vector<Atom> atoms;

  for (int i = 0; i < symbols.size(); i++) {
    Atom a;
    a.symbol = symbols[i];
    a.coord = {coords[3*i], coords[3*i+1], coords[3*i+2]};
    atoms.push_back(a);
  }
  return atoms;
}

std::vector<Fragment> parse_fragments(json input) {
  std::vector<Atom> atoms = parse_atoms(input["molecule"]);

  json fragment_json = input["molecule"]["fragments"];

  int nfrag = fragment_json["nfrag"];

  std::vector<Fragment> fragment_vector(nfrag);
  for (int i = 0; i < nfrag; i++)
    fragment_vector[i].charge = fragment_json["fragment_charges"][i];

  for (int i = 0; i < fragment_json["fragid"].size(); i++) {
    int fragid = fragment_json["fragid"][i];
    fragment_vector[fragid-1].atoms.push_back(atoms[i]);
  }

  return fragment_vector;
}

json encode_topology(std::vector<Fragment> fragments) {
  json topology = json::object();
  topology["geometry"] = json::array();
  topology["symbols"] = json::array();
  topology["fragments"] = json::array();
  topology["fragment_charges"] = json::array();

  topology["connectivity"] = json::array();

  int natoms = 0;
  for (int fragid = 0; fragid < fragments.size(); fragid++) {
    auto &frag = fragments[fragid];

    topology["fragment_charges"].push_back(frag.charge);
    
    std::vector<int> atom_ids;

    for (auto &atom : frag.atoms) {
      atom_ids.push_back(natoms++);

      // Add atom to symbols and geometry
      topology["symbols"].push_back(atom.symbol);
      for (auto &x : atom.coord)
        topology["geometry"].push_back(x);
    }

    topology["fragments"].push_back(atom_ids);
  }
  return topology;
}

int main(int argc, char *argv[]) {
  std::string filename(argv[1]);
  std::cout << "Opening file " << filename << std::endl;
  std::ifstream f(filename);
  json input = json::parse(f);

  json output = template_output();

  std::vector<Fragment> fragments = parse_fragments(input);
  output["topology"] = encode_topology(fragments);

  std::string ofilename = "output-" + filename;
  std::ofstream of(ofilename);
  std::cout << "Writing to file " << ofilename << std::endl;
  of << output.dump(4) << std::endl;
  return 0;
}
