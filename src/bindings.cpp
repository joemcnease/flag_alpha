#include <cstring>

#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>

#include "h2o.hpp"
#include "co2.hpp"
#include "h2o_co2.hpp"
#include "brine.hpp"
#include "gas.hpp"
#include "oil.hpp"
#include "oil_gas.hpp"
#include "heavy_oil.hpp"


namespace nb = nanobind;


NB_MODULE(_flaglib, m) {

    // Modules
    nb::module_ NB_H2O        = m.def_submodule("h2o");
    nb::module_ NB_Brine      = m.def_submodule("brine");
    nb::module_ NB_H2O_CO2    = m.def_submodule("h2o_co2");
    nb::module_ NB_CO2        = m.def_submodule("co2");
    nb::module_ NB_Gas        = m.def_submodule("gas");
    nb::module_ NB_Oil        = m.def_submodule("oil");
    nb::module_ NB_Oil_Gas    = m.def_submodule("oil_gas");
    nb::module_ NB_Oil_CO2    = m.def_submodule("oil_co2");
    nb::module_ NB_Oil_HC_CO2 = m.def_submodule("oil_hc_co2");
    nb::module_ NB_HeavyOil   = m.def_submodule("heavy_oil");

    // H2O module definitions
    NB_H2O.def("velocity", nb::overload_cast<double, double>(&H2O::velocity));
    NB_H2O.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&H2O::velocity));
    // Example of how to return numpy arrays. Not important right now.
    // NB_H2O.def("velocity_npy",
    //     [](const nb::ndarray<double>& P, const nb::ndarray<double>& T) {
    //         std::vector<double> Pvec(P.size());
    //         std::vector<double> Tvec(T.size());
    //         std::copy(P.data(), P.data() + P.size(), Pvec.begin());
    //         std::copy(T.data(), T.data() + T.size(), Tvec.begin());
    //         auto out = H2O::velocity(Pvec, Tvec);
    //         return nb::ndarray<nb::numpy, double>(out.data(), {out.size()}, nb::handle()).cast();
    // }, nb::rv_policy::reference_internal);
    NB_H2O.def("density", nb::overload_cast<double, double>(&H2O::density));
    NB_H2O.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&H2O::density));
    NB_H2O.def("bulk_modulus", nb::overload_cast<double, double>(&H2O::bulk_modulus));
    NB_H2O.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&H2O::bulk_modulus));
    NB_H2O.def("viscosity", nb::overload_cast<double>(&H2O::viscosity));
    NB_H2O.def("viscosity", nb::overload_cast<const std::vector<double>&>(&H2O::viscosity));
    NB_H2O.def("saturated_vapor_pressure", nb::overload_cast<double>(&H2O::saturated_vapor_pressure));
    NB_H2O.def("saturated_vapor_pressure", nb::overload_cast<const std::vector<double>&>(&H2O::saturated_vapor_pressure));
    NB_H2O.def("saturated_vapor_temperature", nb::overload_cast<double>(&H2O::saturated_vapor_temperature));
    NB_H2O.def("saturated_vapor_temperature", nb::overload_cast<const std::vector<double>&>(&H2O::saturated_vapor_temperature));
    NB_H2O.def("gas_solubility", nb::overload_cast<double, double>(&H2O::gas_solubility));
    NB_H2O.def("gas_solubility", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&H2O::gas_solubility));

    // Brine module definitions
    NB_Brine.def("velocity", nb::overload_cast<double, double, double, double, double, double>(&Brine::velocity));
    NB_Brine.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                               const std::vector<double>&, const std::vector<double>&,
                                               const std::vector<double>&, const std::vector<double>&>(&Brine::velocity));
    NB_Brine.def("density", nb::overload_cast<double, double, double, double, double, double>(&Brine::density));
    NB_Brine.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                              const std::vector<double>&, const std::vector<double>&,
                                              const std::vector<double>&, const std::vector<double>&>(&Brine::density));
    NB_Brine.def("bulk_modulus", nb::overload_cast<double, double, double, double, double, double>(&Brine::bulk_modulus));
    NB_Brine.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                   const std::vector<double>&, const std::vector<double>&,
                                                   const std::vector<double>&, const std::vector<double>&>(&Brine::bulk_modulus));
    NB_Brine.def("viscosity", nb::overload_cast<double, double>(&Brine::viscosity));
    NB_Brine.def("viscosity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&Brine::viscosity));
    NB_Brine.def("resistivity", nb::overload_cast<double, double>(&Brine::resistivity));
    NB_Brine.def("resistivity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&Brine::resistivity));
    NB_Brine.def("gas_solubility", nb::overload_cast<double, double, double>(&Brine::gas_solubility));
    NB_Brine.def("gas_solubility", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                     const std::vector<double>&>(&Brine::gas_solubility));


    // H2O+CO2 module definitions
    NB_H2O_CO2.def("velocity", nb::overload_cast<double, double, double>(&H2O_CO2::velocity));
    NB_H2O_CO2.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                 const std::vector<double>&>(&H2O_CO2::velocity));
    NB_H2O_CO2.def("density", nb::overload_cast<double, double, double>(&H2O_CO2::density));
    NB_H2O_CO2.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                const std::vector<double>&>(&H2O_CO2::density));
    NB_H2O_CO2.def("bulk_modulus", nb::overload_cast<double, double, double>(&H2O_CO2::bulk_modulus));
    NB_H2O_CO2.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                     const std::vector<double>&>(&H2O_CO2::bulk_modulus));

    // CO2 module definitions
    NB_CO2.def("pressure", nb::overload_cast<double, double>(&CO2::pressure));
    NB_CO2.def("pressure", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&CO2::pressure));
    NB_CO2.def("velocity", nb::overload_cast<double, double>(&CO2::velocity));
    NB_CO2.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&CO2::velocity));
    NB_CO2.def("density", nb::overload_cast<double, double>(&CO2::density));
    NB_CO2.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&CO2::density));
    NB_CO2.def("bulk_modulus", nb::overload_cast<double, double>(&CO2::bulk_modulus));
    NB_CO2.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&>(&CO2::bulk_modulus));


    // Gas definitions
    nb::module_ NB_EmpiricalModel1999 = NB_Gas.def_submodule("empirical_model_1999");
    nb::module_ NB_GlobalModel        = NB_Gas.def_submodule("global_model");
    nb::module_ NB_LightModel         = NB_Gas.def_submodule("light_model");
    nb::module_ NB_HydrocarbonModel   = NB_Gas.def_submodule("hydrocarbon_model");
    nb::module_ NB_CombinedModel      = NB_Gas.def_submodule("combined_model");

    NB_EmpiricalModel1999.def("velocity", nb::overload_cast<double, double, double>(&Gas::EmpiricalModel1999::velocity));
    NB_EmpiricalModel1999.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                            const std::vector<double>&>(&Gas::EmpiricalModel1999::velocity));

    NB_GlobalModel.def("velocity", nb::overload_cast<double, double, double>(&Gas::GlobalModel::velocity));
    NB_GlobalModel.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                     const std::vector<double>&>(&Gas::GlobalModel::velocity));
    NB_GlobalModel.def("density", nb::overload_cast<double, double, double>(&Gas::GlobalModel::density));
    NB_GlobalModel.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                    const std::vector<double>&>(&Gas::GlobalModel::density));
    NB_GlobalModel.def("bulk_modulus", nb::overload_cast<double, double, double>(&Gas::GlobalModel::bulk_modulus));
    NB_GlobalModel.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                         const std::vector<double>&>(&Gas::GlobalModel::bulk_modulus));
    
    NB_LightModel.def("velocity", nb::overload_cast<double, double, double>(&Gas::LightModel::velocity));
    NB_LightModel.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                    const std::vector<double>&>(&Gas::LightModel::velocity));
    NB_LightModel.def("density", nb::overload_cast<double, double, double>(&Gas::LightModel::density));
    NB_LightModel.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                   const std::vector<double>&>(&Gas::LightModel::density));
    NB_LightModel.def("bulk_modulus", nb::overload_cast<double, double, double>(&Gas::LightModel::bulk_modulus));
    NB_LightModel.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                        const std::vector<double>&>(&Gas::LightModel::bulk_modulus));

    NB_HydrocarbonModel.def("velocity", nb::overload_cast<double, double, double>(&Gas::HydrocarbonModel::velocity));
    NB_HydrocarbonModel.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                          const std::vector<double>&>(&Gas::HydrocarbonModel::velocity));
    NB_HydrocarbonModel.def("density", nb::overload_cast<double, double, double>(&Gas::HydrocarbonModel::density));
    NB_HydrocarbonModel.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                         const std::vector<double>&>(&Gas::HydrocarbonModel::density));
    NB_HydrocarbonModel.def("bulk_modulus", nb::overload_cast<double, double, double>(&Gas::HydrocarbonModel::bulk_modulus));
    NB_HydrocarbonModel.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                              const std::vector<double>&>(&Gas::HydrocarbonModel::bulk_modulus));

    NB_CombinedModel.def("velocity", nb::overload_cast<double, double, double, double, double>(&Gas::CombinedModel::velocity));
    NB_CombinedModel.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                       const std::vector<double>&, const std::vector<double>&,
                                                       const std::vector<double>&>(&Gas::CombinedModel::velocity));
    NB_CombinedModel.def("density", nb::overload_cast<double, double, double, double, double>(&Gas::CombinedModel::density));
    NB_CombinedModel.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                      const std::vector<double>&, const std::vector<double>&,
                                                      const std::vector<double>&>(&Gas::CombinedModel::density));
    NB_CombinedModel.def("bulk_modulus", nb::overload_cast<double, double, double, double, double>(&Gas::CombinedModel::bulk_modulus));
    NB_CombinedModel.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                           const std::vector<double>&, const std::vector<double>&,
                                                           const std::vector<double>&>(&Gas::CombinedModel::bulk_modulus));
    NB_Gas.def("viscosity", nb::overload_cast<double, double, double>(&Gas::viscosity));
    NB_Gas.def("viscosity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                              const std::vector<double>&>(&Gas::viscosity));


    NB_Oil.def("velocity", nb::overload_cast<double, double, double, double, double>(&Oil::velocity));
    NB_Oil.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                             const std::vector<double>&, const std::vector<double>&,
                                             const std::vector<double>&>(&Oil::velocity));
    NB_Oil.def("density", nb::overload_cast<double, double, double, double, double>(&Oil::density));
    NB_Oil.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                            const std::vector<double>&, const std::vector<double>&,
                                            const std::vector<double>&>(&Oil::density));
    NB_Oil.def("bulk_modulus", nb::overload_cast<double, double, double, double, double>(&Oil::bulk_modulus));
    NB_Oil.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                 const std::vector<double>&, const std::vector<double>&,
                                                 const std::vector<double>&>(&Oil::bulk_modulus));
    NB_Oil.def("bubble_point_pressure", nb::overload_cast<double, double, double, double>(&Oil::bubble_point_pressure));
    NB_Oil.def("bubble_point_pressure", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                          const std::vector<double>&,
                                                          const std::vector<double>&>(&Oil::bubble_point_pressure));
    NB_Oil.def("viscosity", nb::overload_cast<double, double, double, double, double>(&Oil::viscosity));
    NB_Oil.def("viscosity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                              const std::vector<double>&, const std::vector<double>&,
                                              const std::vector<double>&>(&Oil::viscosity));


    NB_Oil_Gas.def("velocity", nb::overload_cast<double, double, double, double, double>(&Oil_Gas::velocity));
    NB_Oil_Gas.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                 const std::vector<double>&, const std::vector<double>&,
                                                 const std::vector<double>&>(&Oil_Gas::velocity));
    NB_Oil_Gas.def("density", nb::overload_cast<double, double, double, double, double>(&Oil_Gas::density));
    NB_Oil_Gas.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                const std::vector<double>&, const std::vector<double>&,
                                                const std::vector<double>&>(&Oil_Gas::density));
    NB_Oil_Gas.def("bulk_modulus", nb::overload_cast<double, double, double, double, double>(&Oil_Gas::bulk_modulus));
    NB_Oil_Gas.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
                                                     const std::vector<double>&, const std::vector<double>&,
                                                     const std::vector<double>&>(&Oil_Gas::bulk_modulus));


    // NB_HeavyOil.def("velocity", nb::overload_cast<double, double, double, double, double>(&HeavyOil::velocity));
    // NB_HeavyOil.def("velocity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
    //                                               const std::vector<double>&, const std::vector<double>&,
    //                                               const std::vector<double>&>(&Oil::velocity));
    // NB_HeavyOil.def("density", nb::overload_cast<double, double, double, double, double>(&Oil::density));
    // NB_HeavyOil.def("density", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
    //                                              const std::vector<double>&, const std::vector<double>&,
    //                                              const std::vector<double>&>(&Oil::density));
    // NB_HeavyOil.def("bulk_modulus", nb::overload_cast<double, double, double, double, double>(&Oil::bulk_modulus));
    // NB_HeavyOil.def("bulk_modulus", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
    //                                                   const std::vector<double>&, const std::vector<double>&,
    //                                                   const std::vector<double>&>(&Oil::bulk_modulus));
    // NB_HeavyOil.def("bubble_point_pressure", nb::overload_cast<double, double, double, double>(&Oil::bubble_point_pressure));
    // NB_HeavyOil.def("bubble_point_pressure", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
    //                                                            const std::vector<double>&,
    //                                                            const std::vector<double>&>(&Oil::bubble_point_pressure));
    // NB_HeavyOil.def("viscosity", nb::overload_cast<double, double, double, double, double>(&Oil::viscosity));
    // NB_HeavyOil.def("viscosity", nb::overload_cast<const std::vector<double>&, const std::vector<double>&,
    //                                                const std::vector<double>&, const std::vector<double>&,
    //                                                const std::vector<double>&>(&Oil::viscosity));
}
