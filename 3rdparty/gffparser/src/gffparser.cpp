#include <gffparser.h>
#include <regex>
#include <sstream>
#include <limits>
using namespace gffparser;

double gff_data_t::d_nodata = std::numeric_limits<double>::max();

bool gffparser::utils::check_no_data(const std::string &str)
{
    return str.empty() || str == "." || str == "na" || str == "NA" || str == "N/A" || str == "n/a";
}

std::string gffparser::utils::join_strmap(const std::unordered_map<std::string, std::string> &strmap)
{
    if(strmap.empty()) return "";
    std::string r;
    for(const auto &s: strmap) {
        r += s.first + ":" + s.second + ";";
    }
    r.resize(r.length()-1);
    return r;
}

std::vector<std::string> gffparser::utils::get_fields(const std::string &line, char delim, bool csv_string_format)
{
    std::vector<std::string> r;
    int delim_cnt = 0;
    auto it = line.find(delim);
    size_t prev = 0;
    while(it != std::string::npos) {
        delim_cnt++;
        std::string part = line.substr(prev, it-prev);
        if(csv_string_format && part.length() > 0 && part.at(0) == '"') {
            auto it2 = line.find('"', prev + 1);
            if(it2 != std::string::npos && it2 > it) {
                prev += 1; // "...
                it = it2;
                part = line.substr(prev, it-prev);
                it += 1; // ..."
            }
            else if(part.front() == '"' && part.back() == '"') {
                part = part.substr(1, part.length()-2);
            }
        }
        r.push_back(std::move(part));
        prev = it + 1;
        it = line.find(delim, prev);
    }
    if(prev != line.length())
        r.push_back(line.substr(prev));
    else if(delim_cnt == r.size())
        r.push_back("");
    return r;
}

std::string gffparser::utils::trim(const std::string &str)
{
    auto s = std::find_if_not(str.begin(), str.end(), [](int c){return std::isspace(c);});
    if(s == str.end()) s = str.begin();
    auto e = std::find_if_not(str.rbegin(), str.rend(), [](int c){return std::isspace(c);});
    if(e == str.rend()) e = str.rbegin();
    return std::string(s, e.base());
}

bool gff_attr_t::is_string() const
{
    return std::holds_alternative<std::string>(_value);
}

bool gff_attr_t::is_integer() const
{
    return std::holds_alternative<int64_t>(_value);
}

bool gff_attr_t::is_float() const
{
    return std::holds_alternative<double>(_value);
}

gff_attr_t &gff_attr_t::set_string(const std::string &val)
{
    _value = val;
    return *this;
}

gff_attr_t &gff_attr_t::set_integer(int64_t val)
{
    _value = val;
    return *this;
}

gff_attr_t &gff_attr_t::set_float(double val)
{
    _value = val;
    return *this;
}

const std::string &gff_attr_t::get_string() const noexcept(false)
{
    return std::get<std::string>(_value);
}

const int64_t &gff_attr_t::get_integer() const noexcept(false)
{
    return std::get<int64_t>(_value);
}

const double &gff_attr_t::get_float() const noexcept(false)
{
    return std::get<double>(_value);
}

bool gff_attr_t::empty() const
{
    return _value.index() == std::variant_npos;
}

std::string gff_attr_t::str() const
{
    if(empty()) return "empty";
    auto idx = _value.index();
    switch (idx) {
    case 0:
        return "string=" + get_string();
    case 1:
        return "integer=" + std::to_string(get_integer());
    case 2:
        return "float=" + std::to_string(get_float());
    default:
        break;
    }
    return "unknown type";
}

bool gff_attr_t::eq(const gff_attr_t &v) const
{
    return _value == v._value;
}

bool convert_str_to_uint64(const std::string &ival, uint64_t &oval)
{
    try {
        oval = std::stoll(ival);
        return true;
    }
    catch(...) {}
    return false;
}

bool convert_str_to_double(const std::string &ival, double &oval)
{
    try {
        oval = std::stod(ival);
        return true;
    }
    catch(...) {}
    return false;
}

gff_data_t gff_parser_t::parse_line(const std::string &line, std::string &error, uint64_t linenum, bool onlystrinval)
{
    auto fields = utils::get_fields(line, '\t', false);
    if(fields.size() != static_cast<int>(gff_field_type_t::FIELDSLEN)) {
        error = "line#" + std::to_string(linenum) + " wrong line format (fields count '" +
                 std::to_string(fields.size()) + "' instead of '" +
                 std::to_string(static_cast<int>(gff_field_type_t::FIELDSLEN)) + "')";
        return gff_data_t();
    }

    gff_data_t linedata;
    linedata.linenum = linenum;

    linedata.position.seqid = fields[(int)gff_field_type_t::SEQID];
    linedata.source = fields[(int)gff_field_type_t::SOURCE];
    linedata.type = fields[(int)gff_field_type_t::TYPE];

    std::string start_s = utils::trim(fields[(int)gff_field_type_t::START]);
    std::string end_s = utils::trim(fields[(int)gff_field_type_t::END]);
    std::string score_s = utils::trim(fields[(int)gff_field_type_t::SCORE]);

    if(utils::check_no_data(start_s) || !convert_str_to_uint64(start_s, linedata.position.start)) {
        error = "line#" + std::to_string(linenum) + " wrong 'start' field";
        return gff_data_t();
    }
    if(utils::check_no_data(end_s) || !convert_str_to_uint64(end_s, linedata.position.end)) {
        error = "line#" + std::to_string(linenum) + " wrong 'end' field";
        return gff_data_t();
    }
    if(!utils::check_no_data(score_s) && !convert_str_to_double(score_s, linedata.score)) {
        error = "line#" + std::to_string(linenum) + " wrong 'score' field";
        return gff_data_t();
    }

    std::string strand_s = utils::trim(fields[(int)gff_field_type_t::STRAND]);
    std::string phase_s = utils::trim(fields[(int)gff_field_type_t::PHASE]);

    if(!utils::check_no_data(strand_s)) {
        if(strand_s.size() != 1) {
            error = "line#" + std::to_string(linenum) + " wrong 'strand' field";
            return gff_data_t();
        }
        linedata.strand = strand_s[0];
    }
    if(!utils::check_no_data(phase_s)) {
        if(phase_s.size() != 1) {
            error = "line#" + std::to_string(linenum) + " wrong 'phase' field";
            return gff_data_t();
        }
        linedata.phase = phase_s[0];
    }

    auto attrlist = utils::get_fields(utils::trim(fields[(int)gff_field_type_t::ATTRIBUTES]), ';', false);
    for(const auto &attr: attrlist) {
        auto nmval = utils::get_fields(attr, '=', false);
        if(nmval.size() != 2) {
            error = "line#" + std::to_string(linenum) + " wrong 'attributes' field format";
            return gff_data_t();
        }
        if(onlystrinval) linedata.set_attr(nmval[0], nmval[1]);
        else linedata.set_attr_auto(nmval[0], nmval[1]);
    }
    return linedata;
}

std::string gff_parser_t::dump(const std::string &prefix) const
{
    std::ostringstream ost;
    ost << prefix << "GFF: "
        << std::to_string(_gff_version) << "\n";
    ost << prefix << "FORCE:\n";
    for(const auto &f: _force_types) {
        ost << prefix << "\t" << f.first << ": ";
        switch (f.second) {
        case force_type_t::ft_str:
            ost << "string\n";
            break;
        case force_type_t::ft_int:
            ost << "integer\n";
            break;
        case force_type_t::ft_flt:
            ost << "float\n";
            break;
        }
    }
    ost << prefix << "DATA:\n";
    for(const auto &d: _data)
        ost << prefix << "\t" << d.str() << "\n";
    return ost.str();
}

gff_parser_t &gff_parser_t::operator<<(const std::string &line)
{
    ++_linenum;
    if(line.empty()) return *this; // empty line nothing to do
    if(!_gff_version) {
        std::regex rg(R"(##gff-version\s(\d+))");
        std::smatch sm;
        if(!std::regex_match(line, sm, rg)) {
            _error = "gff-version tag doesn't found";
            return *this;
        }
        _gff_version = std::stoi(sm[1].str());
        if(_gff_version != 3) {
            _error = "incompatible gff-version ('" + sm[1].str() + "' instead of '3')";
            return *this;
        }
    }
    if(line.at(0) == '#') return *this; // this is comment
    if(_thrpool) {
        _thrpool->push(line, _linenum);
        return *this;
    }
    auto linedata = parse_line(line, _error, _linenum, _onlystrval);
    if(!linedata.empty())
        _data.push_back(std::move(linedata));
    return *this;
}

bool gff_parser_t::has_error() const
{
    return !_error.empty();
}

std::string gff_parser_t::error() const
{
    return _error;
}

bool gff_parser_t::empty() const
{
    return _data.empty();
}

uint64_t gff_parser_t::linenum() const
{
    return _linenum;
}

void gff_parser_t::flush()
{
    if(!_data_by_seqid.empty() && (!_thrpool || _thrpool->empty()))
        return;
    if(_thrpool) {
        _thrpool->flush();
        _data.clear();
         auto pooldata = _thrpool->data();
        std::size_t osize = 0;
        for(auto &t: pooldata)
            osize += t.data.size();
        std::sort(pooldata.begin(), pooldata.end(), [](const gff_data_tmp_t &a, const gff_data_tmp_t &b) {return a.linestart < b.linestart;});
        _data.reserve(osize);
        for(auto &t: pooldata)
            _data.insert(_data.end(), std::make_move_iterator(t.data.begin()), std::make_move_iterator(t.data.end()));
        _thrpool->data().clear();
    }
    _data_by_seqid.clear();
    _data_by_source.clear();
    _data_by_type.clear();
    for(auto &d: _data) {
        auto fseqid = _data_by_seqid.find(d.position.seqid);
        auto fsource = _data_by_source.find(d.source);
        auto ftype = _data_by_type.find(d.type);

        if(fseqid == _data_by_seqid.end()) _data_by_seqid[d.position.seqid] = { &d };
        else _data_by_seqid[d.position.seqid].push_back(&d);

        if(fsource == _data_by_source.end()) _data_by_source[d.source] = { &d };
        else _data_by_source[d.source].push_back(&d);

        if(ftype == _data_by_type.end()) _data_by_type[d.type] = { &d };
        else _data_by_type[d.type].push_back(&d);
    }
}

std::size_t gff_parser_t::size() const
{
    return _data.size();
}

void gff_parser_t::setattr_force_str(const std::string &field)
{
    _force_types[field] = force_type_t::ft_str;
}

void gff_parser_t::setattr_force_int(const std::string &field)
{
    _force_types[field] = force_type_t::ft_int;
}

void gff_parser_t::setattr_force_flt(const std::string &field)
{
    _force_types[field] = force_type_t::ft_flt;
}

std::vector<gff_data_t *> gff_parser_t::get_by_pos(const gff_postition_t &position)
{
    if(_data_by_seqid.empty()) flush();
    std::vector<gff_data_t *> r;
    auto it = _data_by_seqid.find(position.seqid);
    if(it == _data_by_seqid.end()) return r;
    if(position.singlepos()) {
        for(auto &data: it->second) {
            if(data->position.in(position.seqid, position.start))
                r.push_back(data);
        }
    }
    else {
        for(auto &data: it->second) {
            if(data->position.intersect(position.seqid, position.start, position.end))
                r.push_back(data);
        }
    }
    return r;
}

std::vector<gff_data_t *> gff_parser_t::get_by_type(const std::string &type)
{
    if(_data_by_type.empty()) flush();
    auto it = _data_by_type.find(type);
    if(it == _data_by_type.end()) return {};
    return it->second;
}

bool check_attrs(const std::vector<gff_attribute_t> &attr, const gff_data_t &data)
{
    bool eq = true;
    for(const auto &a: attr) {
        eq = data.eq_attr(a.name, a.value);
        if(!eq) {
            for(auto &orv: a.orvalues) {
                eq = data.eq_attr(a.name, orv);
                if(eq) break;
            }
        }
        if(!eq) break;
    }
    return eq;
}

std::vector<gff_data_t *> gff_parser_t::get_by_attr(const std::vector<gff_attribute_t> &attr)
{
    std::vector<gff_data_t *> r;
    for(auto &data: _data) {
        if(check_attrs(attr, data))
            r.push_back(&data);
    }
    return r;
}

std::vector<gff_data_t *> gff_parser_t::get_by(const std::string &type, const std::vector<gff_attribute_t> &attr, const gff_postition_t &position)
{
    if(type.empty() && attr.empty() && position.empty()) return {}; // нет параметров
    if(!type.empty() && attr.empty() && position.empty()) return get_by_type(type); // только первый
    if(type.empty() && !attr.empty() && position.empty()) return get_by_attr(attr); // только второй
    if(type.empty() && attr.empty() && !position.empty()) return get_by_pos(position); // только третий
    std::vector<gff_data_t *> r;
    if(!position.empty()) { // есть позиция
        if(_data_by_seqid.empty()) flush();
        auto it = _data_by_seqid.find(position.seqid);
        if(it == _data_by_seqid.end()) return r;
        for(auto &d: it->second) {
            if(position.singlepos()) {
                if(!d->position.in(position.start)) continue;
            }
            else {
                if(!d->position.intersect(position.start, position.end)) continue;
            }
            if(!type.empty() && d->type != type) continue; // позиция + тип
            if(!attr.empty()) { // позиция + тип + аттрибут или позиция + аттрибут
                if(!check_attrs(attr, *d))
                    continue;
            }
            r.push_back(d);
        }
        return r;
    }
    // остался только вариант тип + аттрибут
    if(_data_by_type.empty()) flush();
    auto it = _data_by_type.find(type);
    if(it == _data_by_type.end()) return r;
    for(auto &d: it->second) {
        if(check_attrs(attr, *d)) r.push_back(d);
    }
    return r;
}

std::vector<gff_data_t *> gff_parser_t::get_by(const std::string &type, const std::vector<gff_attribute_t> &attr)
{
    return get_by(type, attr, gff_postition_t());
}

std::vector<gff_data_t *> gff_parser_t::get_by(const std::string &type, const gff_postition_t &position)
{
    return get_by(type, {}, position);
}

std::vector<gff_data_t *> gff_parser_t::get_by(const std::vector<gff_attribute_t> &attr, const gff_postition_t &position)
{
    return get_by("", attr, position);
}

bool gff_data_t::has_attr(const std::string &name) const
{
    return attributes.find(name) != attributes.end();
}

bool gff_data_t::eq_attr(const std::string &name, const gff_attr_t &val) const
{
    auto atri = attributes.find(name);
    if(atri == attributes.end()) return false;
    return atri->second.eq(val);
}

const gff_attr_t &gff_data_t::get_attr(const std::string &name) const
{
    auto fa = attributes.find(name);
    if(fa == attributes.end()) {
        static gff_attr_t empty;
        return empty;
    }
    return fa->second;
}

void gff_data_t::set_attr_auto(const std::string &name, const std::string &val)
{
    bool isdigit = true;
    int dotcnt = 0;
    for(const auto &c: val) {
        if(c == '.') {
            ++dotcnt;
        }
        else if(!std::isdigit(c)) {
            isdigit = false;
            break;
        }
    }
    try {
        if(!isdigit || dotcnt > 1) set_attr(name, val);
        else if(dotcnt == 1) set_attr(name, std::stod(val));
        else set_attr(name, std::stol(val));
        return;
    }
    catch(...) {}
    set_attr(name, val);
}

void gff_data_t::set_attr(const std::string &name, const std::string &val)
{
    attributes[name] = gff_attr_t().set_string(val);
}

void gff_data_t::set_attr(const std::string &name, int64_t val)
{
    attributes[name] = gff_attr_t().set_integer(val);
}

void gff_data_t::set_attr(const std::string &name, double val)
{
    attributes[name] = gff_attr_t().set_float(val);
}

std::string gff_data_t::str() const
{
    std::ostringstream ost;
    ost << "SEQID:" << position.seqid << " SRC:" << source << " TYPE:" << type
        << " START:" << position.start << " END:" << position.end
        << " SCORE:" << (score != d_nodata ? std::to_string(score) : "n/a")
        << " STRAND:" << std::string(&strand, 1)
        << " PHASE:" << std::string(&phase, 1);
    ost << " ATTR:(";
    auto sz = attributes.size();
    for(const auto &a: attributes) {
        --sz;
        ost << a.first << ":" << a.second.str();
        if(sz) ost << ",";
    }
    ost << ")";
    return ost.str();
}

bool gff_data_t::empty() const
{
    return position.empty();
}

bool gff_postition_t::empty() const
{
    return seqid.empty();
}

bool gff_postition_t::singlepos() const
{
    return start == end;
}

bool gff_postition_t::operator==(const gff_postition_t &pos) const
{
    return seqid == pos.seqid && start == pos.start && end == pos.end;
}

gff_postition_t &gff_postition_t::operator=(const gff_postition_t &pos)
{
    seqid = pos.seqid; start = pos.start; end = pos.end;
    return *this;
}

bool gff_postition_t::in(uint64_t pos) const
{
    return start <= pos && end >= pos;
}

bool gff_postition_t::in(const std::string &sid, uint64_t pos) const
{
    return seqid == sid && start <= pos && end >= pos;
}

bool gff_postition_t::intersect(uint64_t spos, uint64_t epos) const
{
    return std::max(spos,start) <= std::min(epos, end);
}

bool gff_postition_t::intersect(const std::string &sid, uint64_t spos, uint64_t epos) const
{
    return seqid == sid && std::max(spos,start) <= std::min(epos, end);
}

void thread_pool_t::thrproc(std::vector<thr_data_t> &&data, std::mutex *mtx,
                            std::vector<gff_data_tmp_t> *out, std::string *error,
                            bool attrval_onlystr)
{
    gff_data_tmp_t outtmp;
    for(auto &d: data) {
        std::string err;
        auto linedata = gff_parser_t::parse_line(d.data, err, d.lnum, attrval_onlystr);
        if(!err.empty()) {
            std::lock_guard<std::mutex> lk(*mtx);
            *error = err;
        }
        if(!linedata.empty())
            outtmp.data.push_back(std::move(linedata));
    }
    if(!outtmp.data.empty()) {
        outtmp.linestart = data.front().lnum;
        outtmp.lineend = data.back().lnum;
        std::lock_guard<std::mutex> lk(*mtx);
        out->push_back(std::move(outtmp));
    }
}

void thread_pool_t::push(const std::string &line, uint64_t linenum)
{
    constexpr std::size_t max_lines_in_pool = 1000;

    std::lock_guard<std::mutex> lk(_tmpdata_mtx);
    _tmpdata.push_back({line, linenum});
    if(_tmpdata.size() > max_lines_in_pool) {
        makeproc(std::move(_tmpdata), false);
        _tmpdata = std::vector<thr_data_t>();
    }
}

void thread_pool_t::flush()
{
    std::lock_guard<std::mutex> lk(_tmpdata_mtx);
    if(!_tmpdata.empty()) {
        makeproc(std::move(_tmpdata), true);
        _tmpdata = std::vector<thr_data_t>();
    }
    for(auto &t: _tpool)
        if(t.joinable()) t.join();
    _tpool.clear();
}

bool thread_pool_t::empty()
{
    std::lock_guard<std::mutex> lk1(_data_mtx);
    std::lock_guard<std::mutex> lk2(_tmpdata_mtx);
    return _tmpdata.empty() && _data.empty() && _tpool.empty();
}

std::vector<gff_data_tmp_t> thread_pool_t::data()
{
    std::vector<gff_data_tmp_t> r;
    std::swap(r, _data);
    return r;
}

void thread_pool_t::makeproc(std::vector<thr_data_t> &&data, bool sync)
{
    if(sync) {
        for(auto &t: _tpool)
            if(t.joinable()) t.join();
        _tpool.clear();
        thread_pool_t::thrproc(std::move(data), &_data_mtx, &_data, _error, _attronlystr);
        return;
    }
    if(_tpool.size() >= _thrmax) {
        if(_tpool.front().joinable())
            _tpool.front().join();
        _tpool.pop_front();
    }
    _tpool.push_back(std::thread(&thread_pool_t::thrproc, std::move(data), &_data_mtx, &_data, _error, _attronlystr));
}

gff_attribute_t &gff_attribute_t::add_value(const std::string &val)
{
    orvalues.push_back(gff_attr_t(val));
    return *this;
}

gff_attribute_t &gff_attribute_t::add_value(int64_t val)
{
    orvalues.push_back(gff_attr_t(val));
    return *this;
}

gff_attribute_t &gff_attribute_t::add_value(double val)
{
    orvalues.push_back(gff_attr_t(val));
    return *this;
}
