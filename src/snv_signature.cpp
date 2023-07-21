/**
 * @file snv_signature.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Implement Single Variation Mutation mutational signature
 * @version 0.1
 * @date 2023-07-21
 * 
 * @copyright Copyright (c) 2023
 * 
 * MIT License
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <array>
#include <cctype>  // toupper
#include <list>
#include <cmath>
#include <limits>

#include "snv_signature.hpp"


namespace Races
{

namespace Passengers
{

namespace SNV
{


char decode_base(const uint8_t& code)
{
    switch(code) {
        case 0:
            return 'C';
        case 1:
            return 'T';
        case 2:
            return 'G';
        case 3:
            return 'A';
        default:
            throw std::domain_error("Unknown code");
    }
}

uint8_t encode_base(const char& base)
{
    switch(base) {
        case 'C':
        case 'c':
            return 0;
        case 'T':
        case 't':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'A':
        case 'a':
            return 3;
        default:
            throw std::domain_error("Unknown base");
    }
}

char read_a_base(std::istream& in)
{
    char symbol;

    in >> symbol;
    switch(symbol) {
        case 'A':
        case 'C':
        case 'G':
        case 'T':
            return symbol;
        case 'a':
        case 'c':
        case 'g':
        case 't':
            return toupper(symbol);
        default:
        {
            std::ostringstream oss;

            oss << "'" << symbol << "' is not a base.";
            throw std::runtime_error(oss.str());
        }
    }
}

std::istream& read_symbol(std::istream& in, const char& symbol)
{
    char read_symbol;

    in >> read_symbol;

    if (read_symbol != symbol) {
        std::ostringstream oss;

        oss << " Expected '" << symbol << "'; got '" << read_symbol <<"'.";
        throw std::runtime_error(oss.str());
    }

    return in;
}

MutationalContext::MutationalContext():
    code(0)
{}

MutationalContext::MutationalContext(const std::string& nucleic_triplet):
    code(0)
{
    if (nucleic_triplet.size()!=3) {
        throw std::domain_error("Expected a nucleic triplet. Got \""
                                + nucleic_triplet +"\"");
    }

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};
    for (unsigned int i=0; i<3; ++i) {
        code = code | (encode_base(nucleic_triplet[i]) << shifts[i]);
    }
}

MutationalContext::MutationalContext(const uint8_t code)
{
    if (code > 63) {
        throw std::domain_error("'"+std::to_string(code)+"' does not correspond"
                                + " to any mutational context.");
    }

    this->code = code;
}

std::string MutationalContext::get_sequence() const
{
    std::string sequence;

    const uint8_t nucleotide_mask = 0x03;
    
    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};
    for (unsigned int i=0; i<3; ++i) {
        auto base = decode_base((code >> shifts[i])&nucleotide_mask);
        sequence.push_back(static_cast<char>(base));
    }

    return sequence;
}

char MutationalContext::get_central_nucleotide() const
{
    return static_cast<char>(decode_base((code >> 4)&0x03));
}

uint8_t MutationalContext::get_complement(const uint8_t& code)
{
    const uint8_t nucleotide_mask = 0x03;

    // this is to have the central nucleotide associated
    // to the most significant bits
    const std::array<uint8_t, 3> shifts{0,4,2};

    uint8_t complementary_code{0};
    for (unsigned int i=0; i<3; ++i) {
        auto base_code = (code >> shifts[i])&nucleotide_mask;

        // Complement base code by using this tricky feature
        // of base encoding
        base_code =  (base_code + 2) % 4;

        complementary_code = complementary_code | (base_code << shifts[i]);
    }    

    return complementary_code;
}

bool MutationalContext::is_a_base(const char& base)
{
    switch(base) {
        case 'A':
        case 'a':
        case 'C':
        case 'c':
        case 'G':
        case 'g':
        case 'T':
        case 't':
            return true;
        default:
            return false;
    }
}

char MutationalContext::get_complement(const char& base)
{
    switch(base) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';
        default:
        {
            std::ostringstream oss;

            oss << "Unsupported base '" << base << "'";
            throw std::domain_error(oss.str());
        }
    }
}

std::string MutationalContext::get_complement(const std::string& sequence)
{
    std::string complementary(sequence);

    for (char& nucleotide: complementary) {
        nucleotide = get_complement(nucleotide);
    }

    return complementary;
}

MutationalType::MutationalType():
    context(0), replace_base('A')
{}

MutationalType::MutationalType(const MutationalContext& context, const char& replace_base)
{
    auto central_nucleotide = context.get_central_nucleotide();

    if (central_nucleotide == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context.get_sequence() +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!MutationalContext::is_a_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (central_nucleotide != 'C' && central_nucleotide != 'T') {
        this->context = context.get_complement(); 
        this->replace_base = MutationalContext::get_complement(replace_base);
    } else {
        this->context = context;
        this->replace_base = toupper(replace_base);
    }
}

MutationalType::MutationalType(const std::string& context, const char& replace_base):
    context(context), replace_base(toupper(replace_base))
{
    if (context[1] == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!MutationalContext::is_a_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (context[1] != 'C' &&  context[1] != 'T' && context[1] != 'c' &&  context[1] != 't') {
        this->context = this->context.get_complement(); 
        this->replace_base = MutationalContext::get_complement(replace_base);
    }
}

MutationalType::MutationalType(const std::string& type)
{
    std::istringstream iss(type);

    iss >> *this;
}

MutationalSignatureExprValue::MutationalSignatureExprValue(const std::map<MutationalType, double>& dist_map):
    value_map(dist_map)
{}

MutationalSignatureExprValue::operator MutationalSignature()
{
    return MutationalSignature(value_map);
}

MutationalSignatureExprValue& MutationalSignatureExprValue::operator+(MutationalSignatureExprValue&& expression_value)
{
    for (auto& [type, prob]: value_map) {
        prob += expression_value.value_map[type];
    }

    return *this;
}

MutationalSignatureExprValue& MutationalSignatureExprValue::operator+(const MutationalSignature& signature)
{
    for (auto& [type, prob]: value_map) {
        prob += signature(type);
    }

    return *this;
}

MutationalSignature::MutationalSignature():
    dist_map()
{
    dist_map[MutationalType()] = 1;
}

bool is_about(double x, double y)
{
    return std::fabs(x - y) <= std::numeric_limits<double>::epsilon();
}

MutationalSignature::MutationalSignature(const std::map<MutationalType, double>& distribution):
    dist_map(distribution)
{
    double partial = 0;
    for (const auto& [type, prob]: dist_map) {
        if (prob<0 || prob>1) {
            std::ostringstream oss;

            oss << type << ": " << prob;
            throw std::domain_error("The parameter is not a probability distribution: "
                                    + oss.str());
        }
        partial += prob;
    }
    if (!is_about(partial,1)) {
        std::ostringstream oss;
        oss << "The parameter is not a probability distribution: the sum of"
            << " the probabilities is " << partial;

        throw std::domain_error(oss.str());
    }
}

std::vector<std::string> read_row(std::istream& in, const char& delimiter)
{
    std::string line;
    getline(in, line);

    std::istringstream iss(line);
    
    std::vector<std::string> row;
    std::string cell;

    while (getline(iss, cell, delimiter)) {
        row.push_back(cell);
    }

    return row;
}

std::map<std::string, std::map<MutationalType, double>> not_validated_read_from_stream(std::istream& in, const char& delimiter)
{
    std::map<std::string, std::map<MutationalType, double>> result;
    std::vector<std::string> name_vector = read_row(in,delimiter);

    if (name_vector.front() != "") {
        std::ostringstream oss;

        oss << "The cell in position (1,1) should be an empty string. Read \"" << name_vector.front() << "\".";
        throw std::runtime_error(oss.str());
    }

    unsigned int row_number = 2;
    std::vector<std::string> row = read_row(in,delimiter);
    while (row.size() != 0) {
        if (row.size() != name_vector.size()) {

            std::cout << row.size() << std::endl;
            std::cout << name_vector.size() << std::endl;
            std::ostringstream oss;

            oss << "The header and the row number " << row_number << " differ in size.";
            throw std::runtime_error(oss.str());
        }

        MutationalType type(row.front());

        auto row_it = row.begin();
        auto name_it = name_vector.begin();
        for (++row_it, ++name_it;row_it != row.end(); ++row_it, ++name_it) {
            result[*name_it][type] = std::stod(*row_it);
        }

        ++row_number;
        row = read_row(in,delimiter);
    } 

    return result;
}


std::map<std::string, MutationalSignature> MutationalSignature::read_from_stream(std::istream& in)
{
    std::map<std::string, MutationalSignature> result;

    auto name_signature_map = not_validated_read_from_stream(in,'\t');
    for (const auto& [name, signature]: name_signature_map) {
        try {
            result.emplace(std::string(name), MutationalSignature(signature));
        } catch (std::domain_error& ex) {
            std::ostringstream oss;

            oss << "Column \"" << name << "\" is not a mutational signature: "
                << ex.what();
            throw std::runtime_error(oss.str());
        }
    }

    return result;
}

}  // SNV

}  // Passengers

}  // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Passengers::SNV::MutationalType& type)
{
    std::string type_sequence = type.get_context().get_sequence();

    out << type_sequence[0] << "[" 
        << type_sequence[1] << ">" << type.get_replace_base() 
        << "]" << type_sequence[2];

    return out;
}

std::istream& operator>>(std::istream& in, Races::Passengers::SNV::MutationalType& type)
{
    using namespace Races::Passengers::SNV;

    std::string seq;
    char replace_base;

    seq.push_back(read_a_base(in));
    read_symbol(in, '[');
    seq.push_back(read_a_base(in));
    read_symbol(in, '>');
    in >> replace_base;
    read_symbol(in, ']');
    seq.push_back(read_a_base(in));

    try {
        type = MutationalType(seq, replace_base);
    } catch (std::domain_error& err) {
        throw std::runtime_error(err.what());
    }
    return in;
}

}   // std
