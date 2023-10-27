/**
 * @file driver.cpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Implements the driver genotype representation
 * @version 0.10
 * @date 2023-10-28
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

#include <iostream>
#include <sstream>

#include "driver_genotype.hpp"
#include "cell_event.hpp"

namespace Races 
{

namespace Drivers
{

unsigned int EpigeneticGenotype::counter = 0;
unsigned int Genotype::counter = 0;

EpigeneticRates::EpigeneticRates(const double methylation_rate, const double demethylation_rate):
        methylation(methylation_rate), demethylation(demethylation_rate)
{
    if (methylation_rate<0 || methylation_rate>1) {
        throw std::domain_error("the methylation rate does not belong to the interval [0,1]");
    }

    if (demethylation_rate<0 || demethylation_rate>1) {
        throw std::domain_error("the demethylation rate does not belong to the interval [0,1]");
    }
}

EpigeneticRates::EpigeneticRates(const double rate):
        EpigeneticRates(rate, rate)
{}

EpigeneticRates& EpigeneticRates::set_methylation_rate(const double& rate)
{
    if (rate<0 || rate>1) {
        throw std::domain_error("the rate does not belong to the interval [0,1]");
    }

    methylation = rate;

    return *this;
}

EpigeneticRates& EpigeneticRates::set_demethylation_rate(const double& rate)
{
    if (rate<0 || rate>1) {
        throw std::domain_error("the rate does not belong to the interval [0,1]");
    }

    demethylation = rate;

    return *this;
}

EpigeneticGenotype::EpigeneticGenotype():
    id(0), genomic_id(0)
{}

EpigeneticGenotype::EpigeneticGenotype(const Genotype& genotype,
                                       const size_t num_of_promoters):
    id(counter++), genomic_id(genotype.get_id())
{
    const size_t epigenetic_index = genotype.epigenetic_genotypes().size();

    name = genotype.get_name();

    methylation_signature = Genotype::index_to_signature(epigenetic_index, num_of_promoters);
}

double EpigeneticGenotype::get_rate(const CellEventType& event) const
{
    if (event_rates.count(event)==0) {
        return 0;
    }
    return event_rates.at(event);
}

std::string EpigeneticGenotype::get_epigenetic_name() const
{
    return name + Genotype::signature_to_string(methylation_signature);
}

Genotype::Genotype(const std::string& name,
                   const std::vector<EpigeneticRates>& epigenetic_event_rates):
    id(counter++), name(name)
{
    size_t epigenetic_mutations = 1<<epigenetic_event_rates.size();

    for (size_t i=0; i<epigenetic_mutations; ++i) {
        e_genotypes.push_back(EpigeneticGenotype(*this, epigenetic_event_rates.size()));
    }

    for (auto& e_genotype: e_genotypes) {
        auto e_signature = e_genotype.get_methylation_signature();
        auto& e_rates = e_genotype.epigenetic_rates;
        
        for (size_t i=0; i<epigenetic_event_rates.size(); ++i) {
            // get the signature of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            e_signature[i] = !e_signature[i];

            // get the index of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            size_t index = signature_to_index(e_signature);

            // get the identifier of the genotype reachable by 
            // methylating/demethylating the i-th promoter
            EpigeneticGenotypeId dst_id = e_genotypes[index].get_id();

            if (e_signature[i]) {  // if the promoter is methylated
                // a methylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].get_methylation_rate();
            } else {
                // a demethylation event must occur
                e_rates[dst_id] = epigenetic_event_rates[i].get_demethylation_rate();
            }

            // revert to the original signature
            e_signature[i] = !e_signature[i];
        }
    }
}

Genotype::Genotype(const std::string& name):
    Genotype(name, std::vector<EpigeneticRates>())
{}

size_t Genotype::num_of_promoters() const
{
    if (e_genotypes.size()==0) {
        return 0;
    }

    return e_genotypes[0].get_methylation_signature().size();
}

EpigeneticGenotype& Genotype::operator[](const MethylationSignature& methylation_signature)
{
    validate_signature(methylation_signature);

    return e_genotypes[Genotype::signature_to_index(methylation_signature)];
}

const EpigeneticGenotype& Genotype::operator[](const MethylationSignature& methylation_signature) const
{
    validate_signature(methylation_signature);

    return e_genotypes[Genotype::signature_to_index(methylation_signature)];
}

EpigeneticGenotype& Genotype::operator[](const std::string& methylation_signature)
{
    validate_signature(methylation_signature);

    return e_genotypes[Genotype::string_to_index(methylation_signature)];
}

const EpigeneticGenotype& Genotype::operator[](const std::string& methylation_signature) const
{
    validate_signature(methylation_signature);

    return e_genotypes[Genotype::string_to_index(methylation_signature)];
}

size_t Genotype::string_to_index(const std::string& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const char& c: methylation_signature) {
        digit_value >>= 1;
        if (c=='+') {
            index |= digit_value;
        }
    }

    return index;
}

std::string Genotype::index_to_string(const size_t& index, const size_t num_of_promoters)
{
    size_t digit_value=(1 << num_of_promoters);
    std::ostringstream oss;

    while (digit_value>1) {
        digit_value >>= 1;
        oss << (((digit_value&index) != 0) ? "+" : "-");
    }

    return oss.str();
}

size_t Genotype::signature_to_index(const MethylationSignature& methylation_signature)
{
    size_t index=0, digit_value=(1 << methylation_signature.size());

    for (const bool& b_value: methylation_signature) {
        digit_value >>= 1;
        if (b_value) {
            index |= digit_value;
        }
    }

    return index;
}

MethylationSignature Genotype::index_to_signature(const size_t& index, const size_t num_of_promoters)
{
    MethylationSignature signature;
    size_t digit_value=(1 << num_of_promoters);

    while (digit_value>1) {
        digit_value >>= 1;
        signature.push_back((digit_value&index) != 0);
    }

    return signature;
}

}   // Drivers

}   // Races

namespace std
{

std::ostream& operator<<(std::ostream& out, const Races::Drivers::EpigeneticRates& epigentic_rates)
{
    out << "{\"on\": " << epigentic_rates.get_methylation_rate() 
            << ",\"off\": " << epigentic_rates.get_demethylation_rate() << "}";
    return out;
}

std::ostream& operator<<(std::ostream& out, const Races::Drivers::EpigeneticGenotype& genotype)
{
    out << "{name: \""<< genotype.get_epigenetic_name() << "\", id: " << genotype.get_id() 
        << ", event_rates: {";

    std::string sep="";
    for (const auto& [event, rate]: genotype.get_rates()) {
        out << sep << Races::Drivers::cell_event_names[event] << ": " << rate;
        sep = ", ";
    }

    const auto& e_rates = genotype.get_epigenetic_mutation_rates();
    if (e_rates.size()>0) {
        out << "}, epigenetic rates: {";
        sep="";
        for (const auto& [dst_id, rate]: e_rates) {
            out << sep << dst_id << ": " << rate;
            sep = ", ";
        }
        out << "}";
    }

    out << "}";

    return out;
}

std::ostream& operator<<(std::ostream& out, const Races::Drivers::Genotype& genotype)
{
    out << "{name: \"" << genotype.get_name() << "\", id: " 
        << genotype.get_id() << ", epigenetic genotypes=[";

    const auto& e_genotypes = genotype.epigenetic_genotypes();
    
    std::string sep = "";
    if (e_genotypes.size()>1) {
        sep = "\n";
    }

    for (const auto& e_genotype: e_genotypes) {
        out << sep << e_genotype;
        sep = ",\n";
    }

    out << "]}";
    return out; 
}

}  // std