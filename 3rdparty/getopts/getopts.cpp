#include "getopts.h"
#include <filesystem>
#include <sstream>

get_opts_t::get_opts_t(int argc, char **argv)
{
    auto argv_path = std::filesystem::path(argv[0]);
    _program_path = std::filesystem::absolute(argv_path).string();
    _program_name = argv_path.filename().string();
    cmd_opt_t *curopt = nullptr;
    for(int i = 1; i < argc; ++i) {
        std::string v(argv[i]);
        if(v.empty()) continue;
        if(v.size() > 1 && v.at(0) == '-') { // parameter
            auto optname = v.substr(1);
            curopt = nullptr;
            for(auto &o: _opts) {
                if(o.name == optname) {
                    curopt = &o;
                    break;
                }
            }
            if(!curopt) {
                _opts.push_back(cmd_opt_t());
                curopt = &_opts.back();
                curopt->name = optname;
            }
        }
        else {
            if(!curopt) {
                _opts.push_back(cmd_opt_t());
                curopt = &_opts.back();
            }
            curopt->values.push_back(v);
        }
    }
}

std::vector<cmd_opt_t> &get_opts_t::result()
{
    return _opts;
}

std::string &get_opts_t::program_name()
{
    return _program_name;
}

std::string &get_opts_t::program_path()
{
    return _program_path;
}

std::string get_opts_t::dump(const std::string &prefix) const
{
    std::ostringstream ost(std::ios::ate);
    for(const auto &o: _opts) {
        ost << prefix << o.name << "\n";
        for(const auto &v: o.values)
            ost << prefix << "  " << v << "\n";
    }
    return ost.str();
}

bool cmd_opt_t::equal(const std::string &optname, bool canstartpart) const
{
    if(optname.size()) {
        if(canstartpart && name.find(optname) == 0)
            return true;
        return name == optname;
    }
    return name.empty();
}
