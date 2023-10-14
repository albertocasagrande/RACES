/**
 * @file json_config.hpp
 * @author Alberto Casagrande (alberto.casagrande@uniud.it)
 * @brief Defines classes and function for reading JSON configurations
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

#ifndef __RACES_JSON_CONFIG__
#define __RACES_JSON_CONFIG__

#include <map>
#include <vector>

#include <nlohmann/json.hpp>

#include "simulation.hpp"
#include "position_set.hpp"
#include "mutation_engine.hpp"

namespace Races
{

namespace Passengers 
{

class ConfigReader
{
    /**
     * @brief Get a sample region corner
     * 
     * @param sampler_region_json is the JSON of the region to sample
     * @param corner_field_name is the corner field name
     * @return std::vector<Races::Drivers::Simulation::AxisPosition> 
     */
    static std::vector<Races::Drivers::Simulation::AxisPosition>
    get_corner(const nlohmann::json& sampler_region_json, const std::string& corner_field_name);

    /**
     * @brief Get a mutation rate
     * 
     * @param mutation_rate_json is the JSON of the mutation rate
     * @return a pair epigenetic status/mutation rate
     */
    static std::pair<std::string, double>
    get_mutation_rate(const nlohmann::json& mutation_rate_json);
    
    /**
     * @brief Get the CNA type
     * 
     * @param CNA_json is the JSON of the CNA
     * @return the type of the CNA
     */
    static Races::Passengers::CopyNumberAlteration::Type
    get_CNA_type(const nlohmann::json& CNA_json);

    /**
     * @brief Add a CNA to a list
     * 
     * @param CNAs is the list of CNA in which the read CNA should be inserted
     * @param CNA_json is the JSON of the CNA to be inserted into the list
     */
    static void
    add_CNA(std::list<Races::Passengers::CopyNumberAlteration>& CNAs, 
            const nlohmann::json& CNA_json);

    /**
     * @brief Add a SNV to a list
     * 
     * @param SNVs is the list of SNV in which the read SNV should be inserted
     * @param SNV_json is the JSON of the SNV to be inserted into the list
     */
    static void
    add_SNV(const std::string& name, std::list<Races::Passengers::SNV>& SNVs,
            const nlohmann::json& SNV_json);

    /**
     * @brief Get the mutation rates
     * 
     * @param mutation_rates_json is the JSON of the mutation rates
     * @return a map associating an epigenetic status to its mutation rate
     */
    static std::map<std::string, double>
    get_mutation_rates(const nlohmann::json& mutation_rates_json);

    /**
     * @brief Add the driver mutations to SNV and CNA lists
     * 
     * @param driver_name is the driver name
     * @param SNVs is the list of SNVs in which the new SNVs must be inserted
     * @param CNAs is the list of CNAs in which the new CNAs must be inserted
     * @param driver_mutation_json is the JSON of the driver mutation
     */
    static void
    add_driver_mutation(const std::string& driver_name,
                        std::list<Races::Passengers::SNV>& SNVs,
                        std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                        const nlohmann::json& driver_mutation_json);

    /**
     * @brief Get the fraction field
     * 
     * @param fraction_json is the JSON of the faction field
     * @return the fraction value 
     */
    static double get_fraction(const nlohmann::json& fraction_json);

    /**
     * @brief Add driver mutational properties
     * 
     * @param mutational_properties is the species mutational properties
     * @param drivers_simulation is the driver simulation
     * @param driver_properties_json is the JSON of the driver properties
     */
    static void 
    add_driver_mutational_properties(SpeciesMutationalProperties& mutational_properties,
                                     const Races::Drivers::Simulation::Simulation& drivers_simulation,
                                     const nlohmann::json& driver_properties_json);

public:
    /**
     * @brief Get the default mutational coefficients
     * 
     * @param mutational_coeff_json is the JSON of the mutational coefficients
     * @return a map associating to a set of mutational signatures their percentage
     *          in the default mutational configuration 
     */
    static std::map<std::string, double>
    get_default_mutational_coefficients(const nlohmann::json& mutational_coeff_json);

    /**
     * @brief Add the timed mutational coefficients to a mutation engine 
     * 
     * @tparam GENOME_WIDE_POSITION is the genome wide position type
     * @tparam RANDOM_GENERATOR is the type of the random generator
     * @param engine is a mutation engine
     * @param mutational_coeff_json is the JSON of the mutational coefficients
     */
    template<typename GENOME_WIDE_POSITION, typename RANDOM_GENERATOR>
    static void
    add_timed_mutational_coefficients(Races::Passengers::MutationEngine<GENOME_WIDE_POSITION,RANDOM_GENERATOR>& engine, 
                                      const nlohmann::json& mutational_coeff_json)
    {
        std::map<double, std::map<std::string, double>> timed_coefficients;

        if (!mutational_coeff_json.is_object()) {
            throw std::runtime_error("The \"mutational coefficients\" field must be an object");
        }

        if (mutational_coeff_json.contains("timed")) {

            auto& timed_coeff_json = mutational_coeff_json["timed"];

            if (!timed_coeff_json.is_array()) {
                throw std::runtime_error("The optional \"timed\" field must be an array");
            }

            for (const auto& timed_json : timed_coeff_json) {
                if (!timed_json.is_object()) {
                    throw std::runtime_error("The elements of the \"timed\" field must be objects");
                }

                double time = get_from<double>("time", timed_json, 
                                               "The elements of the \"timed\" field");

                if (timed_coefficients.count(time)>0) {
                    throw std::runtime_error("Two elements of the \"timed\" field have the "
                                             "same \"time\"");
                }

                expecting("coefficients", timed_json, "The elements of the \"timed\" field");

                engine.add(time, get_mutational_coefficients(mutational_coeff_json["default"]));
            }
        }
    }

    /**
     * @brief Get the sample set
     * 
     * @param sample_region_json is the sample region JSON
     * @return the rectangle position set
     */
    static Drivers::RectangleSet 
    get_sample_region(const nlohmann::json& sample_region_json);

    /**
     * @brief Collect the driver mutations
     * 
     * @param[in] driver_name is the driver name
     * @param[in,out] SNVs is the list of driver SNVs
     * @param[in,out] CNAs is the list of driver CNAs
     * @param[in] driver_mutations_json is the JSON of driver mutations
     */
    static void
    collect_driver_mutations(const std::string& driver_name,
                             std::list<Races::Passengers::SNV>& SNVs,
                             std::list<Races::Passengers::CopyNumberAlteration>& CNAs,
                             const nlohmann::json& driver_mutations_json);

    /**
     * @brief Get the species mutational properties
     * 
     * @param drivers_simulation is a driver simulation
     * @param simulation_json is the JSON of the simulation configuration
     * @return the species mutational properties
     */
    static Races::Passengers::SpeciesMutationalProperties
    get_mutational_properties(const Races::Drivers::Simulation::Simulation& drivers_simulation,
                              const nlohmann::json& simulation_json);

    /**
     * @brief Raise an exception if a field is not available in a JSON
     * 
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     */
    static 
    void expecting(const std::string& field_name, const nlohmann::json& json, 
                   const std::string& context="");

    /**
     * @brief Get a value of a JSON field
     * 
     * @tparam TYPE is the type of the JSON field value
     * @param field_name is the aimed field name
     * @param json is the JSON that must contain the field
     * @param context is a description of the context
     * @return the value of the JSON field `field_name`
     */
    template<typename TYPE>
    static TYPE get_from(const std::string& field_name, const nlohmann::json& json, 
                         const std::string& context="")
    {
        expecting(field_name, json, context);

        return json[field_name].template get<TYPE>();
    } 

    /**
     * @brief Get the number of wild-type genome alleles from the configuration
     * 
     * @param simulation_json is the JSON simulation configuration
     * @return the number of wild-type genome alleles
     */
    static uint8_t
    get_number_of_alleles(const nlohmann::json& simulation_json);

    /**
     * @brief Get the mutational coefficients
     * 
     * @param mutational_coefficients_json the JSON of the mutational coefficients
     * @return a map associating a set of mutational signature to their percentage
     *          in the mutational coefficients
     */
    static std::map<std::string, double>
    get_mutational_coefficients(const nlohmann::json& mutational_coefficients_json);
};

}   // Drivers

}   // Races

#endif // __RACES_JSON_CONFIG__