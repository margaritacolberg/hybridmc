#include "main_helpers.h"
namespace po = boost::program_options;

// wang landau plugin for python function
std::vector<double> wang_landau_process(const std::string& json_name, std::optional<std::string> input_name) {

    std::ifstream input(json_name);
    nlohmann::json json;
    input >> json;

    const Param p = json;
    const Box box{p.length};
    UpdateConfig update_config;

    const unsigned int t_bonds = p.transient_bonds.get_nbonds();
    const unsigned int nbonds = p.nonlocal_bonds.get_nbonds();
    const unsigned int nstates = std::pow(2, t_bonds); // NOLINT(cppcoreguidelines-narrowing-conversions)
    ConfigInt store_config_int;
    std::vector<uint64_t> config_count(nstates);
    std::vector<double> dist(nbonds);

    // initialize both members of CountBond struct to 0
    CountBond count_bond = {};

    System sys(p.nbeads);

    std::seed_seq seq(p.seeds.begin(), p.seeds.end());
    Random mt(seq);

    if (input_name.has_value() && std::filesystem::exists(*input_name)) {
        read_input(*input_name, sys.pos);
        init_update_config(sys.pos, update_config, box, p.transient_bonds);
        init_s(sys.s_bias, t_bonds);
    } else {
        init_pos(sys.pos, box, mt, p);
        init_s(sys.s_bias, t_bonds);
    }

    if (!check_local_dist(sys.pos, box, p.near_min2, p.near_max2, p.nnear_min2,
                          p.nnear_max2)) {
        throw std::runtime_error("local beads overlap");
    }

    if (!check_nonlocal_dist(sys.pos, box, p.rh2, p.stair2,
                             p.transient_bonds, p.permanent_bonds)) {
        throw std::runtime_error("nonlocal beads overlap");
    }

    //get bead index for transient pair and rc value for it
    std::tuple<unsigned int, unsigned int, double> transient_pair = p.transient_bonds.getBond(0);
    int bead1 = std::get<0>(transient_pair);
    int bead2 = std::get<1>(transient_pair);
   // std::cout << "bead1: " << bead1 << " bead2: " << bead2 << std::endl;
   // the minimum rc2 can be is the target bonding rc for this nonlocal bond; get this value and assign to rc_min2
    double rc_min2 = std::get<2>(transient_pair);
    // take square root to get the minimum bonding distance value
    double rc_min = std::sqrt(rc_min2);
    // distance vector to store all distances
    std::vector<double> distance_values;

    double gamma = p.gamma;
    unsigned int iter_wl = 0;
    unsigned int native_ind = nstates - 1;
    double wall_time = 0.0;

    while (gamma > p.gamma_f) {
        // iterate over the 2 states
        for (unsigned int i = 0; i < nstates; i++) {
            // run trajectory to get final state
            Config state = run_trajectory_wl(sys, mt, p, box, update_config,
                                             count_bond, wall_time, iter_wl);

            if (state == 0) {
                sys.s_bias[native_ind] -= gamma;
            } else {
                sys.s_bias[native_ind] += gamma;
            }

            double dx = sys.pos[bead2].x - sys.pos[bead1].x;
            double dy = sys.pos[bead2].y - sys.pos[bead1].y;
            double dz = sys.pos[bead2].z - sys.pos[bead1].z;
            box.mindist(dx, dy, dz);

            const double distance = sqrt(dx * dx + dy * dy + dz * dz);

            //std::cout << " recording everything" << dist << " sb = " << sys.s_bias[native_ind] << std::endl;

            if (distance > rc_min){
                //std::cout << " recording " << dist << " sb = " << sys.s_bias[native_ind] << std::endl;
                distance_values.push_back(distance);
            }
        }

        iter_wl += 1;
        gamma = 1.0 / double(iter_wl);
    }

    //sort distance vector
    std::sort(distance_values.begin(), distance_values.end());

    std::cout << "sbias is " << sys.s_bias[0] - sys.s_bias[1] << std::endl;
    std::cout << "native sbias is " << sys.s_bias[1] << std::endl;
    std::cout << "non native sbias is " << sys.s_bias[0] << std::endl;

    // create return vector
    std::vector<double> return_info;
    // get s bias value
    return_info.push_back(sys.s_bias[0] - sys.s_bias[1]);
    // get rc_min
    return_info.push_back(rc_min);
    // get size of distance value array
    int n_rcs = distance_values.size();
    // want integer index of n_rcs of 0.10 (10%), 0.30 (30%) etc
    double rc_middle = distance_values[int(n_rcs * 0.005)];
    double rc_outermost = distance_values[int(n_rcs * 0.05 )];

    std::cout << "n_rcs: " << n_rcs << " rc_middle: " << rc_middle << " rc_outermost: " << rc_outermost << std::endl;
    std::cout << "True rc " << rc_min << " rc_min: " << distance_values[0] << " rc_max: " << distance_values[n_rcs - 1] << std::endl;

    return_info.push_back(rc_middle);
    return_info.push_back(rc_outermost);

  return return_info;
}

// create pybind11 module for wang_landau function
namespace py = pybind11;
PYBIND11_MODULE(wang_landau, m) {
    m.doc() = "pybind11 hybridmc plugin for wang_landau function";

    using namespace py::literals;
    m.def("WL_process", &wang_landau_process, "A function that conducts the wang landau algorithm",
          "json_name"_a, "h5_input_name"_a=py::none());

}