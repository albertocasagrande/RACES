/**
 * @file json_config.cpp
 * @author Alberto Casagrande (alberto.casagrande@units.it)
 * @brief Implements classes and function for reading JSON configurations
 * @version 0.2
 * @date 2023-10-14
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

#include "json_config.hpp"

namespace Races
{

namespace Passengers 
{

std::vector<Races::Drivers::Simulation::AxisPosition>
ConfigReader::get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name)
{
    if (!sampler_region_json.contains(corner_field_name)) {
        throw std::runtime_error("The \"sampler region\" must contain a "
                                    "\"" + corner_field_name + "\" field");
    }

    auto corner_json = sampler_region_json[corner_field_name];

    if (!corner_json.is_array()) {
        throw std::runtime_error("The \"" + corner_field_name + "\" must be an array of natural values");
    }

    std::vector<Races::Drivers::Simulation::AxisPosition> corner;

    for (const auto& axis_value : corner_json) {
        corner.push_back(axis_value.template get<Races::Drivers::Simulation::AxisPosition>());
    }

    if (corner.size()<2 || corner.size()>3) {
        throw std::runtime_error("\""+corner_field_name+"\" must be either a 2D or 3D array");
    }

    return corner;
}

std::pair<std::string, double>
ConfigReader::get_mutation_rate(const nlohmann::json& mutation_rate_json)
{  
    if (!mutation_rate_json.is_object()) {
        throw std::runtime_error("All the elements in \"mutation rates\" "
                                    "must be objects");
    }

    return {
        get_from<std::string>("epigenetic status", mutation_rate_json, 
                                "All the elements in \"mutation rates\""),
        get_from<double>("rate", mutation_rate_json, 
                            "All the elements in \"mutation rates\"")
    };
}

Races::Passengers::CopyNumberAlteration::Type
ConfigReader::get_CNA_type(const nlohmann::json& CNA_json)
{
    auto subtype_str = get_from<std::string>("subtype", CNA_json, "All the CNAs");
    for (auto& subtype_chr: subtype_str) {
        subtype_chr = tolower(subtype_chr);
    }

    if (subtype_str == "amplification") {
        return Races::Passengers::CopyNumberAlteration::Type::AMPLIFICATION;
    }

    if (subtype_str == "deletion") {
        return Races::Passengers::CopyNumberAlteration::Type::DELETION;
    }

    throw std::runtime_error("Unknown CNA \""+
                                CNA_json["subtype"].template get<std::string>()+"\"");
}

void
ConfigReader::add_CNA(std::list<Races::Passengers::CopyNumberAlteration>& CNAs, 
                      const nlohmann::json& CNA_json)
{
    using namespace Races::Passengers;

    auto chr_str = get_from<std::string>("chromosome", CNA_json, "All the CNAs");
    auto position = get_from<ChrPosition>("position", CNA_json, "All the CNAs");
    auto length = get_from<ChrPosition>("length", CNA_json, "All the CNAs");
    auto allele = get_from<uint8_t>("allele", CNA_json, "All the CNAs");

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    switch (get_CNA_type(CNA_json)) {
        case CopyNumberAlteration::Type::AMPLIFICATION:
            CNAs.push_back(CopyNumberAlteration::new_amplification({genomic_position, length}, 
                                                                   allele, 0));
            break;
        case CopyNumberAlteration::Type::DELETION:
            CNAs.push_back(CopyNumberAlteration::new_deletion({genomic_position, length}, allele));
            break;
        default:
            throw std::runtime_error("Unsupported CNA type");
    }
}

void
ConfigReader::add_SNV(const std::string& driver_name, std::list<Races::Passengers::SNV>& SNVs,
                      const nlohmann::json& SNV_json)
{
    using namespace Races::Passengers;

    MutationalContext context = get_from<std::string>("context", SNV_json, 
                                                        "All the SNVs");
    auto chr_str = get_from<std::string>("chromosome", SNV_json, "All the SNVs");
    auto mutated_base = get_from<std::string>("mutated base", SNV_json, "All the SNVs");
    auto position = get_from<ChrPosition>("position", SNV_json, "All the SNVs");

    GenomicPosition genomic_position(GenomicPosition::stochr(chr_str), position);

    SNVs.emplace_back(genomic_position, context, mutated_base[0], driver_name);
}

std::map<std::string, double>
ConfigReader::get_mutation_rates(const nlohmann::json& mutation_rates_json)
{
    std::map<std::string, double> mutational_rates;

    if (!mutation_rates_json.is_array()) {
        throw std::runtime_error("The \"mutation rates\" field must be an array "
                                "of driver mutation rates");
    }

    for (const auto& mutation_rate_json : mutation_rates_json) {
        mutational_rates.insert(get_mutation_rate(mutation_rate_json));
    }

    return mutational_rates;
}

void
ConfigReader::add_driver_mutation(const std::string& driver_name,
                                  std::list<Races::Passengers::SNV>& SNVs,
                                  std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                                  const nlohmann::json& driver_mutation_json)
{
    if (!driver_mutation_json.is_object()) {
        throw std::runtime_error("All the elements in \"mutations\" must be objects");
    }

    auto type = get_from<std::string>("type", driver_mutation_json,
                                        "All the elements in \"mutations\"");
    if (type=="SNV") {
        add_SNV(driver_name, SNVs, driver_mutation_json);

        return;
    }

    if (type=="CNA") {
        add_CNA(CNAs, driver_mutation_json);

        return;
    }

    throw std::runtime_error("Unsupported mutations type \""+type+"\"");
}

double
ConfigReader::get_fraction(const nlohmann::json& fraction_json)
{
    const double fraction = fraction_json.template get<double>();

    if (fraction<0) {
        throw std::runtime_error("The \"fraction\" field must be a non-negative value");
    }

    return fraction;
}

void
ConfigReader::add_driver_mutational_properties(Races::Passengers::SpeciesMutationalProperties& mutational_properties,
                                               const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                               const nlohmann::json& driver_properties_json)
{
    if (!driver_properties_json.is_object()) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                    "must be objects");
    }

    auto driver_name = get_from<std::string>("name", driver_properties_json, 
                                             "All the elements in \"driver properties\"");

    if (!driver_properties_json.contains("mutation rates")) {
        throw std::runtime_error("All the elements in \"driver properties\" "
                                "must contain a \"mutation rates\" field");
    }
    auto mutation_rates = get_mutation_rates(driver_properties_json["mutation rates"]);

    using namespace Races::Passengers;
    std::list<SNV> SNVs;
    std::list<CopyNumberAlteration> CNAs;
    if (driver_properties_json.contains("mutations")) {
        collect_driver_mutations(driver_name, SNVs, CNAs, driver_properties_json["mutations"]);
    }
    
    mutational_properties.add_species(drivers_simulation, driver_name, mutation_rates, SNVs, CNAs);
}

std::map<std::string, double>
ConfigReader::get_default_mutational_coefficients(const nlohmann::json& mutational_coeff_json)
{
    if (!mutational_coeff_json.is_object()) {
        throw std::runtime_error("The \"mutational coefficients\" field must be an object");
    }

    expecting("default", mutational_coeff_json, "The \"mutational coefficients\" field");

    return get_mutational_coefficients(mutational_coeff_json["default"]);
}

Drivers::RectangleSet
ConfigReader::get_sample_region(const nlohmann::json& sampler_region_json)
{
    if (!sampler_region_json.is_object()) {
        throw std::runtime_error("The \"sample region\" field must be objects");
    }

    auto lower_vector = get_corner(sampler_region_json, "lower corner");
    auto upper_vector = get_corner(sampler_region_json, "upper corner");

    if (lower_vector.size() != upper_vector.size()) {
        throw std::runtime_error("The \"lower corner\" and \"upper corner\" arrays must "
                                    "have the same size");
    }

    using namespace Races::Drivers::Simulation;

    if (lower_vector.size()==2) {
        return {{lower_vector[0], lower_vector[1]},
                {upper_vector[0], upper_vector[1]}};
    }

    return {{lower_vector[0], lower_vector[1], lower_vector[2]},
            {upper_vector[0], upper_vector[1], upper_vector[2]}};
}

void 
ConfigReader::collect_driver_mutations(const std::string& driver_name,
                                        std::list<Races::Passengers::SNV>& SNVs,
                                        std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                                        const nlohmann::json& driver_mutations_json)
{
    if (!driver_mutations_json.is_array()) {
        throw std::runtime_error("The \"mutations\" field must be an array "
                                    "of mutations");
    }

    for (const auto& driver_mutation_json : driver_mutations_json) {
        add_driver_mutation(driver_name, SNVs, CNAs, driver_mutation_json);
    }
}

Races::Passengers::SpeciesMutationalProperties
ConfigReader::get_mutational_properties(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                        const nlohmann::json& simulation_json)
{
    using namespace Races::Passengers;

    SpeciesMutationalProperties mutational_properties;

    expecting("driver properties", simulation_json, "The passengers simulation configuration");

    auto& mutational_properties_json = simulation_json["driver properties"];

    if (!mutational_properties_json.is_array()) {
        throw std::runtime_error("The \"driver properties\" field must be an array "
                                    "of driver mutational properties");
    }

    for (const auto& driver_properties_json : mutational_properties_json) {
        add_driver_mutational_properties(mutational_properties, drivers_simulation,
                                         driver_properties_json);
    }

    return mutational_properties;
}


void ConfigReader::expecting(const std::string& field_name, const nlohmann::json& json, 
                             const std::string& context)
{
    if (!json.contains(field_name)) {
        if (context.length()>0) {
            throw std::runtime_error("Expected a missing \""+field_name+"\" field");
        } else {
            throw std::runtime_error(context+ " must contains a \""+field_name+"\" field");
        }
    }
}

uint8_t
ConfigReader::get_number_of_alleles(const nlohmann::json& simulation_json)
{
    const double value = get_from<uint8_t>("number of alleles", simulation_json,
                                           "The passengers simulation configuration");

    if (value <= 0) {
        throw std::runtime_error("The number of allele must be positive");
    }

    return value;
}

std::map<std::string, double>
ConfigReader::get_mutational_coefficients(const nlohmann::json& mutational_coefficients_json)
{
    if (!mutational_coefficients_json.is_array()) {
        throw std::runtime_error("Mutational coefficients must be an array of object");
    }

    double total=0;

    std::map<std::string, double> mutational_coefficients;
    for (const auto& mutational_coefficient_json : mutational_coefficients_json) {

        if (!mutational_coefficient_json.is_object()) {
            throw std::runtime_error("Mutational coefficient must be an object");
        }

        const std::string SBS = get_from<std::string>("SBS", mutational_coefficient_json,
                                                        "Every mutational coefficient");
        if (mutational_coefficients.count(SBS)>0) {
            throw std::runtime_error("\""+SBS+"\" already among mutational coefficients");
        }

        double fraction = get_from<double>("fraction", mutational_coefficient_json,
                                            "Every mutational coefficient");

        total += fraction;
        if (total>1) {
            std::ostringstream oss;

            oss << "The sum of the \"fraction\" fields must be 1. "
                << "The difference is " << (1-total);

            throw std::runtime_error(oss.str());
        }

        mutational_coefficients[SBS] = fraction;
    }
    
    if (abs(1-total)>10*std::numeric_limits<double>::epsilon()) {
        std::ostringstream oss;

        oss << "The sum of the \"fraction\" fields must be 1. "
            << "The difference is " << (1-total);

        throw std::runtime_error(oss.str());
    }

    return mutational_coefficients;
}

}   // Passengers

}   // Races
