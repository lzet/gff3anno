#include <filesystem>
#include <iostream>
#include <gffparser.h>
#include <getopts.h>
#include <bxzstr.hpp>
#include <regex>

enum class finput_type_t {
    fi_err = -2, fi_unk = -1, fi_bed = 0, fi_vcf = 1
};

void usage(const std::string &program) {
    std::cerr << "USAGE: "
              << "\n         gff column types: seqid, source, type, pos, endpos, score, strand, phase, attr\n"
              << "\n         " << program
              << "\n         -gff [path/to/file.gff3] #input gff3 file"
              << "\n         -type {bed,vcf} #input file type (default: get from extension)"
              << "\n         -skip N #(optional) for bed-file only, skip lines (default 0)"
              << "\n         -header N #for bed-file only, header line num (default -1, no header)"
              << "\n         -threads N #max threads (default 4)"
              << "\n         -in path/to/input.{bed,vcf} #input bed or vcf file (or '-' for stdin, '-type' required)"
              << "\n         -out path/to/output.{bed,vcf} #output file path or '-' for stdout"
              << "\n         -seqid N #sequence id column number (default 1)"
              << "\n         -pos N #position column number (default 2)"
              << "\n         -endpos N #(optional) end position column number (info:<name> for vcf)"
              << "\n         -where <par1>...<parN> #(optional) select from gff parameter (format <coltype>[:<attrname>]:<value>)"
              << "\n         -add <par1>...<parN> #fields to add to output file (format: <coltype>[:<attrname>])"
              << "\n         -ext {intersect,length} #add extended information ('intersect' - intersect percent)"
              << "\n"
              << std::endl;
    std::cerr << "USAGE EXAMPLE: "
              << "\n         " << program << " -h"
              << "\n         " << program << " -gff gencode.v47.primary_assembly.basic.annotation.gff3 -in test1.bed -out - -seqid 1 -pos 2 -where type:gene attr:gene_name:ADA -add attr:gene_name attr:gene_id type"
              << "\n"
              << std::endl;
}

int arg_error(const std::string &parname, const std::string &program) {
    usage(program);
    std::cerr << "'" << parname << "' parameter must be set" << std::endl;
    return 1;
}

int noarg_error(const std::string &parname, const std::string &program) {
    usage(program);
    std::cerr << "parameter '" << parname << "' is not yet implemented" << std::endl;
    return 1;
}

int arg_type_error(const std::string &parname, const std::string &error, const std::string &program) {
    usage(program);
    std::cerr << "error in '" << parname << "': " << error << std::endl;
    return 1;
}

int get_colnum(const std::vector<std::string> &v) {
    if(v.size() == 1) {
        try {
            return std::stoi(v.front());
        }
        catch (...) {}
    }
    return -1;
}

enum class select_column_t {
    noval = -1,
    seqid = 0,
    source,
    type,
    pos,
    endpos,
    score,
    strand,
    phase,
    attr,

    intersect,
    length
};

struct selectpar_t {
    select_column_t colnum = select_column_t::noval;
    std::string attrname;
    std::string value;
    static select_column_t colnum_by_name(const std::string &msg) {
        if(msg == "seqid") return select_column_t::seqid;
        if(msg == "source") return select_column_t::source;
        if(msg == "type") return select_column_t::type;
        if(msg == "pos") return select_column_t::pos;
        if(msg == "endpos") return select_column_t::endpos;
        if(msg == "score") return select_column_t::score;
        if(msg == "strand") return select_column_t::strand;
        if(msg == "phase") return select_column_t::phase;
        if(msg == "attr") return select_column_t::attr;
        if(msg == "intersect") return select_column_t::intersect;
        if(msg == "length") return select_column_t::length;
        return select_column_t::noval;
    }
    static std::string colname_by_num(select_column_t type) {
        switch (type) {
        case select_column_t::noval:  return "noval";
        case select_column_t::seqid:  return "seqid";
        case select_column_t::source: return "source";
        case select_column_t::type:   return "type";
        case select_column_t::pos:    return "pos";
        case select_column_t::endpos: return "endpos";
        case select_column_t::score:  return "score";
        case select_column_t::strand: return "strand";
        case select_column_t::phase:  return "phase";
        case select_column_t::attr:   return "attr";
        case select_column_t::intersect:   return "intersect";
        case select_column_t::length:   return "length";
        }
        return "noval";
    }
    std::string init() {
        auto fields = gffparser::utils::get_fields(_format, ':', false);
        if(_novalue) {
            if(fields.size() == 1) {
                colnum = colnum_by_name(fields[0]);
            }
            else if(fields.size() == 2) {
                colnum = colnum_by_name(fields[0]);
                attrname = fields[1];
            }
            if(colnum == select_column_t::noval)
                return "unknown format for '" + _format + "' (for example: 'type' or 'attr:gene_name')";
        }
        else {
            if(fields.size() == 2) {
                colnum = colnum_by_name(fields[0]);
                value = fields[1];
            }
            else if(fields.size() == 3) {
                colnum = colnum_by_name(fields[0]);
                attrname = fields[1];
                value = fields[2];
            }
            if(colnum == select_column_t::noval)
                return "unknown format for '" + _format + "' (for example: 'type:gene' or 'attr:gene_name:ADA')";
        }
        return "";
    }
    std::string orig() const {
        return _format;
    }
    selectpar_t(const std::string &txt, bool novalue)
        : _format(txt), _novalue(novalue) {}
private:
    std::string _format;
    bool _novalue;
};

std::string join_strlist(const std::vector<std::string> &strlist, char delim = ',')
{
    if(strlist.empty()) return "";
    std::string r;
    for(const auto &s: strlist) {
        r += s + delim;
    }
    r.resize(r.length()-1);
    return r;
}

finput_type_t check_extension(const std::string &filepath, bool &gzipped)
{
    std::regex fext_r(R"(.*(\.bed|\.vcf)(\.gz)?$)");
    std::smatch fm;
    if(std::regex_match(filepath, fm, fext_r)) {
        if(fm.size() >= 2) {
            auto ext = fm[1].str();
            gzipped = fm.size() > 2;
            if(ext == ".bed") return finput_type_t::fi_bed;
            if(ext == ".vcf") return finput_type_t::fi_vcf;
        }
    }
    return finput_type_t::fi_err;
}

std::string intersect_percent(uint64_t s1, uint64_t e1, uint64_t s2, uint64_t e2)
{
    int64_t ival = std::min(e1, e2) - std::max(s1, s2) + 1;
    if(ival <= 0) return "0";
    int64_t l1 = e1 - s1 + 1;
    if(l1 == 0) throw std::runtime_error("gel length error (=0)");
    int64_t p = (ival * 1000) / l1;
    return std::to_string(p/10) + "." + std::to_string(p%10);
}

int main(int argc, char **argv)
{
    get_opts_t inopts(argc, argv);
    std::filesystem::path gffpath;
    bxz::ifstream ifile;
    std::ofstream ofile;
    std::istream *iptr = nullptr;
    std::ostream *optr = nullptr;
    int seqid = 1;
    int pos = 2;
    int endpos = -1;
    std::string endpos_vcf;
    int nproc = 4;
    int skip = 0;
    int header = -1;
    bool gzipped = false;
    finput_type_t ftype = finput_type_t::fi_unk;
    std::vector<selectpar_t> where, add;
    std::string ifpath, ofpath;
    for(auto &p: inopts.result()) {
        if(p.equal("h")) {
            usage(inopts.program_name());
            return 0;
        }
        if(p.equal("gff")) {
            if(p.values.size() == 1)
                gffpath = p.values.front();
            continue;
        }
        if(p.equal("in")) {
            if(p.values.size() == 1) {
                ifpath = p.values.front();
                if(ifpath == "-")
                    iptr = &std::cin;
                else if(ifpath != ofpath)
                    ifile.open(ifpath);
                if(!iptr && ifile.is_open())
                    iptr = &ifile;
            }
            continue;
        }
        if(p.equal("out")) {
            if(p.values.size() == 1) {
                ofpath = p.values.front();
                if(ofpath == "-")
                    optr = &std::cout;
                else if(ifpath != ofpath)
                    ofile.open(ofpath);
                if(!optr && ofile.is_open())
                    optr = &ofile;
            }
            continue;
        }
        if(p.equal("threads")) {
            nproc = get_colnum(p.values);
            continue;
        }
        if(p.equal("seqid")) {
            seqid = get_colnum(p.values);
            continue;
        }
        if(p.equal("pos")) {
            pos = get_colnum(p.values);
            continue;
        }
        if(p.equal("endpos")) {
            endpos = get_colnum(p.values);
            if(endpos < 0 && p.values.size() == 1) {
                if(p.values.front().size() > 5 && p.values.front().substr(0, 5) == "info:")
                    endpos_vcf = p.values.front().substr(5);
            }
            continue;
        }
        if(p.equal("where")) {
            for(const auto &v: p.values)
                where.push_back(selectpar_t(v, false));
            continue;
        }
        if(p.equal("add")) {
            for(const auto &v: p.values)
                add.push_back(selectpar_t(v, true));
            continue;
        }
        if(p.equal("skip")) {
            skip = get_colnum(p.values);
            continue;
        }
        if(p.equal("header")) {
            header = get_colnum(p.values);
            continue;
        }
        if(p.equal("type")) {
            if(p.values.size() == 1 && p.values.front() == "bed")
                ftype = finput_type_t::fi_bed;
            else if(p.values.size() == 1 && p.values.front() == "vcf")
                ftype = finput_type_t::fi_vcf;
            else
                ftype = finput_type_t::fi_err;
            continue;
        }
        if(p.equal("ext")) {
            for(const auto &v: p.values) {
                if(v == "intersect") {
                    add.push_back(selectpar_t("intersect", true));
                    continue;
                }
                if(v == "length") {
                    add.push_back(selectpar_t("length", true));
                    continue;
                }
            }
            continue;
        }
    }
    if(ftype == finput_type_t::fi_unk)
        ftype = check_extension(ifpath, gzipped);
    else check_extension(ifpath, gzipped);

    if(gffpath.empty() || !std::filesystem::exists(gffpath))
        return arg_error("-gff", inopts.program_name());
    if(ofpath == ifpath)
        return arg_error("-in/-out same", inopts.program_name());
    if(!iptr)
        return arg_error("-in", inopts.program_name());
    if(!optr)
        return arg_error("-out", inopts.program_name());
    if(nproc < 0)
        return arg_error("-threads", inopts.program_name());
    if(seqid < 0)
        return arg_error("-seqid", inopts.program_name());
    if(pos < 0)
        return arg_error("-pos", inopts.program_name());
    if(ftype == finput_type_t::fi_err || ftype == finput_type_t::fi_unk)
        return arg_error("-type", inopts.program_name());
    for(auto &w: where) {
        auto er = w.init();
        if(!er.empty())
            return arg_type_error("-where", er, inopts.program_name());
    }
    for(auto &a: add) {
        auto er = a.init();
        if(!er.empty())
            return arg_type_error("-add", er, inopts.program_name());
    }

    gffparser::gff_parser_t gff(nproc, true);
    {
        bxz::ifstream gffifs(gffpath);
        if(!gffifs.is_open()) {
            usage(inopts.program_name());
            std::cerr << "can't open '" << gffpath.string() << "'" << std::endl;
            return 1;
        }
        std::string line;
        while(std::getline(gffifs, line)) {
            gff << line;
            if(gff.has_error()) {
                std::cerr << "[GFF ERROR] " << gff.error() << std::endl;
                return 1;
            }
        }
        gff.flush();
        if(gff.has_error()) {
            std::cerr << "[GFF ERROR] " << gff.error() << std::endl;
            return 1;
        }
    }
    std::string inln;
    if(endpos < 0) // single pos mod
        endpos = pos;
    seqid -= 1; pos -= 1; endpos -= 1;
    std::vector<gffparser::gff_attribute_t> attr_v;
    std::string type_s;
    bool hasval = false;
    for(const auto &w: where) {
        switch (w.colnum) {
        case select_column_t::attr:
            hasval = false;
            for(auto &a: attr_v) {
                if(a.name == w.attrname) {
                    hasval = true;
                    a.orvalues.push_back(gffparser::gff_attr_t(w.value));
                    break;
                }
            }
            if(!hasval) attr_v.push_back(gffparser::gff_attribute_t(w.attrname, w.value));
            break;
        case select_column_t::type:
            if(type_s.empty()) type_s = w.value;
            else return arg_error("-type must be only one", inopts.program_name());
            break;
        case select_column_t::noval:
        case select_column_t::seqid:
        case select_column_t::source:
        case select_column_t::pos:
        case select_column_t::endpos:
        case select_column_t::score:
        case select_column_t::strand:
        case select_column_t::phase:
            return noarg_error("-where " + w.orig(), inopts.program_name());
        case select_column_t::intersect:
        case select_column_t::length:
            // ext field, nothing to do
            break;
        }
    }

    std::string *seqid_s = nullptr;
    uint64_t pos_ui = 0;
    uint64_t endpos_ui = 0;

    auto getansw = [&]() -> std::vector< std::vector<std::string> > {
        auto items = gff.get_by(type_s, attr_v, gffparser::gff_postition_t(*seqid_s, pos_ui, endpos_ui));
        std::vector< std::vector<std::string> > addfields;
        addfields.resize(add.size());
        for(const auto &i: items) {
            for(std::size_t ai = 0; ai < add.size(); ++ai) {
                auto &a = add[ai];
                auto &to = addfields[ai];
                switch (a.colnum) {
                case select_column_t::noval: break;
                case select_column_t::seqid: to.push_back(i->position.seqid); break;
                case select_column_t::source: to.push_back(i->source); break;
                case select_column_t::type: to.push_back(i->type); break;
                case select_column_t::pos: to.push_back(std::to_string(i->position.start)); break;
                case select_column_t::endpos: to.push_back(std::to_string(i->position.end)); break;
                case select_column_t::score: to.push_back(std::to_string(i->score)); break;
                case select_column_t::strand: to.push_back(std::to_string(i->strand)); break;
                case select_column_t::phase: to.push_back(std::to_string(i->phase)); break;
                case select_column_t::attr: to.push_back(i->get_attr(a.attrname).get_string()); break;
                case select_column_t::intersect: to.push_back(intersect_percent(pos_ui, endpos_ui, i->position.start, i->position.end)); break;
                case select_column_t::length: to.push_back(std::to_string(endpos_ui-pos_ui+1)); break;
                }
            }
        }
        return addfields;
    };

    try {
        std::size_t lncnt = 0;
        if(ftype == finput_type_t::fi_bed) {
            while(std::getline(*iptr, inln)) {
                ++lncnt;
                std::ostringstream ost;
                ost << inln;
                if(header > 0 && lncnt == header) {
                    header = 0;
                    for(const auto &a: add)
                        ost << "\t" << a.orig();
                    ost << "\n";
                }
                else if(skip > 0) {
                    --skip;
                    ost << "\n";
                }
                else if(inln.empty() || inln[0] == '#') {
                    ost << "\n";
                }
                else {
                    auto bedfields = gffparser::utils::get_fields(inln, '\t', false);
                    seqid_s = &bedfields.at(seqid);
                    pos_ui = std::stol(bedfields.at(pos));
                    endpos_ui = pos != endpos ? std::stol(bedfields.at(endpos)) : pos;
                    for(const auto &af: getansw()) {
                        ost << "\t" << join_strlist(af);
                    }
                    ost << "\n";
                }
                std::string ost_s = ost.str();
                optr->write(ost_s.c_str(), ost_s.size());
            }
        }
        else if(ftype == finput_type_t::fi_vcf) {
            bool need_info = true;
            while(std::getline(*iptr, inln)) {
                ++lncnt;
                std::ostringstream ost;
                if(inln.empty() || inln[0] == '#') {
                    ost << inln << "\n";
                    if(need_info && inln.size() > 7 && inln.substr(0, 7) == "##INFO=") {
                        need_info = false;
                        for(const auto &a: add)
                            ost << "##INFO=<ID=" << a.attrname
                                << ",Number=1,Type=String,Description=\""
                                << inopts.program_name() << " " << a.attrname << "\">\n";
                    }
                }
                else {
                    auto tabi1 = inln.find('\t');
                    auto tabi0 = 0;
                    // #CHROM-0  POS-1 ID-2  REF-3 ALT-4 QUAL-5    FILTER-6  INFO-7    FORMAT-8  sample-name-9
                    uint32_t intpos = 0;
                    std::string seqid_stmp;
                    seqid_s = &seqid_stmp;
                    pos_ui = 0;
                    endpos_ui = 0;
                    while(tabi1 != std::string::npos) {
                        if(intpos == 0) { // chrom
                            seqid_stmp = inln.substr(0, tabi1);
                        }
                        else if(intpos == 1) {
                            std::string pos_s = inln.substr(tabi0 + 1, tabi1 - tabi0 - 1);
                            pos_ui = std::stol(pos_s);
                        }
                        else if(intpos == 7) { // info
                            if(endpos_vcf.empty()) {
                                endpos_ui = pos_ui;
                            }
                            else {
                                auto ei = inln.find(endpos_vcf, tabi0 + 1); // name=<val>;
                                if(ei != std::string::npos) {
                                    ei += endpos_vcf.size() + 1; // sizeof(name=)
                                    auto endpos_s = inln.substr(ei, inln.find_first_of(";\t", ei + 1) - ei);
                                    endpos_ui = std::stol(endpos_s);
                                }
                                else {
                                    endpos_ui = pos_ui;
                                }
                            }
                            if(tabi1 > 2 && inln[tabi1-1] == '.' && inln[tabi1-2] == '\t') { // INFO = .
                                ost << inln.substr(0, tabi1-1);
                                auto af = getansw();
                                for(std::size_t i = 0; i < add.size(); ++i) {
                                    if(i) ost << ";";
                                    ost << add[i].attrname << "=" << join_strlist(af[i]);
                                }
                                ost << inln.substr(tabi1);
                            }
                            else {
                                ost << inln.substr(0, tabi1) << ";";
                                auto af = getansw();
                                for(std::size_t i = 0; i < add.size(); ++i) {
                                    if(i) ost << ";";
                                    ost << add[i].attrname << "=" << join_strlist(af[i]);
                                }
                                ost << inln.substr(tabi1);
                            }
                        }
                        tabi0 = tabi1;
                        tabi1 = inln.find('\t', tabi1+1);
                        ++intpos;
                    }
                    ost << "\n";
                }
                std::string ost_s = ost.str();
                optr->write(ost_s.c_str(), ost_s.size());
            }
        }
    }
    catch(const std::exception &e) {
        std::cerr << "[" << inopts.program_name() << " ERROR] " << e.what() << "\n";
        std::cerr << "columns: seqid=" << (seqid+1) << " pos=" << (pos+1) << " endpos=" << (endpos+1) << "\n";
        std::cerr << "last value: seqid=" << (seqid_s ? *seqid_s : "null") << " pos=" << (pos_ui) << " endpos=" << (endpos_ui) << "\n";
        std::cerr << "last line: '" << inln << "'" << std::endl;
    }

    return 0;
}
