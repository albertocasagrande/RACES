/**
 * @file sbs_signature.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements SBS signature
 * @version 0.17
 * @date 2024-05-11
 *
 * @copyright Copyright (c) 2023-2024
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

#include "sbs_signature.hpp"


namespace Races
{

namespace Mutations
{

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

SBSType::SBSType():
    context(), replace_base('A')
{}

SBSType::SBSType(const SBSContext& context, const char& replace_base)
{
    auto central_nucleotide = context.get_central_nucleotide();

    if (central_nucleotide == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context.get_sequence() +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!GenomicSequence::is_a_DNA_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (central_nucleotide != 'C' && central_nucleotide != 'T') {
        this->context = context.get_complemented();
        this->replace_base = GenomicSequence::get_complemented(replace_base);
    } else {
        this->context = context;
        this->replace_base = toupper(replace_base);
    }
}

SBSType::SBSType(const std::string& context, const char& replace_base):
    context(context), replace_base(toupper(replace_base))
{
    if (context[1] == replace_base) {
        std::ostringstream oss;

        oss << "Expected a replace base different from "
            << " the second nucleotide in the context. Got \""
            << context +"\" and '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (!GenomicSequence::is_a_DNA_base(replace_base)) {
        std::ostringstream oss;

        oss << "Expected a replace base. Got '" << replace_base << "'";

        throw std::domain_error(oss.str());
    }

    if (context[1] != 'C' &&  context[1] != 'T' && context[1] != 'c' &&  context[1] != 't') {
        this->context = this->context.get_complemented();
        this->replace_base = GenomicSequence::get_complemented(replace_base);
    }
}

SBSType::SBSType(const std::string& type)
{
    std::istringstream iss(type);

    iss >> *this;
}

SBSSignatureExprResult::SBSSignatureExprResult():
    value_map()
{}

SBSSignatureExprResult::SBSSignatureExprResult(const std::map<SBSType, double>& dist_map):
    value_map(dist_map)
{}

SBSSignatureExprResult::operator SBSSignature()
{
    return SBSSignature(value_map);
}

SBSSignatureExprResult& SBSSignatureExprResult::operator+(SBSSignatureExprResult&& expression_value)
{
    for (auto& [type, prob]: value_map) {
        auto it = expression_value.value_map.find(type);
        if (it != expression_value.value_map.end()) {
            prob += expression_value.value_map[type];
        }
    }

    for (const auto& [type, prob]: expression_value.value_map) {
        if (value_map.count(type)==0) {
            value_map[type] = prob;
        }
    }

    return *this;
}

SBSSignatureExprResult& SBSSignatureExprResult::operator+(const SBSSignature& signature)
{
    for (auto& [type, prob]: value_map) {
        prob += signature(type);
    }

    for (const auto& [type, prob]: signature) {
        if (value_map.count(type)==0) {
            value_map[type] = prob;
        }
    }

    return *this;
}

SBSSignature::SBSSignature():
    dist_map()
{
    dist_map[SBSType()] = 1;
}

bool is_about(double x, double y)
{
    return std::fabs(x - y) <= 1e-6;
}

SBSSignature::SBSSignature(const std::map<SBSType, double>& distribution):
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
        oss << "The parameter is not a probability distribution: 1 minus the sum of"
            << " the probabilities is " << (1-partial);

        throw std::domain_error(oss.str());
    }
}

double SBSSignature::operator()(const SBSType& type) const
{
    auto it = dist_map.find(type);
    if (it != dist_map.end()) {
        return it->second;
    }

    return 0;
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

std::map<std::string, std::map<SBSType, double>> not_validated_read_from_stream(std::istream& in, const char& delimiter)
{
    std::map<std::string, std::map<SBSType, double>> result;
    std::vector<std::string> name_vector = read_row(in,delimiter);

    unsigned int row_number = 2;
    std::vector<std::string> row = read_row(in,delimiter);
    while (row.size() != 0) {
        if (row.size() != name_vector.size()) {
            std::ostringstream oss;

            oss << "The header and the row number " << row_number << " differ in size.";
            throw std::runtime_error(oss.str());
        }

        SBSType type(row.front());

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


std::map<std::string, SBSSignature> SBSSignature::read_from_stream(std::istream& in)
{
    std::map<std::string, SBSSignature> result;

    std::string curLocale = setlocale(LC_ALL, nullptr);
    setlocale(LC_ALL,"C");
    auto name_signature_map = not_validated_read_from_stream(in,'\t');
    setlocale(LC_ALL, curLocale.c_str());
    for (const auto& [name, signature]: name_signature_map) {
        try {
            result.emplace(std::string(name), SBSSignature(signature));
        } catch (std::domain_error& ex) {
            std::ostringstream oss;

            oss << "Column \"" << name << "\" is not a SBS signature: "
                << ex.what();
            throw std::runtime_error(oss.str());
        }
    }

    return result;
}

std::map<std::string, SBSSignature> SBSSignature::read_from_stream(std::istream& in, const std::set<std::string>& signature_names)
{
    auto result = SBSSignature::read_from_stream(in);

    auto it = result.begin();

    while (it != result.end()) {
        if (signature_names.count(it->first)==0) {
            result.erase(it++);
        } else {
            ++it;
        }
    }

    return result;
}

}  // Mutations

}  // Races

namespace std
{

bool less<Races::Mutations::SBSType>::operator()(const Races::Mutations::SBSType &lhs,
                                                         const Races::Mutations::SBSType &rhs) const
{
    const auto& lhs_code = lhs.get_context().get_code();
    const auto& rhs_code = rhs.get_context().get_code();

    return ((lhs_code < rhs_code) ||
            ((lhs_code == rhs_code) && (lhs.get_replace_base()<rhs.get_replace_base())));
}

std::ostream& operator<<(std::ostream& out, const Races::Mutations::SBSType& type)
{
    std::string type_sequence = type.get_context().get_sequence();

    out << type_sequence[0] << "["
        << type_sequence[1] << ">" << type.get_replace_base()
        << "]" << type_sequence[2];

    return out;
}

std::istream& operator>>(std::istream& in, Races::Mutations::SBSType& type)
{
    using namespace Races::Mutations;

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
        type = SBSType(seq, replace_base);
    } catch (std::domain_error& err) {
        throw std::runtime_error(err.what());
    }
    return in;
}

}   // std
