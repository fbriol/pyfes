#include "fes/settings.hpp"

#include <boost/range/algorithm/find.hpp>
#include <vector>
#include <map>

#include "fes/constituent.hpp"
#include "fes/interface/wave.hpp"
#include "fes/interface/wave_table.hpp"

namespace fes {

inline auto check_if_is_inferred(
    const ConstituentId ident, const bool is_modeled,
    const std::vector<ConstituentId>& inferred_constituents) -> bool {
  if (is_modeled) {
    return false;
  }
  return boost::range::find(inferred_constituents, ident) !=
         inferred_constituents.end();
}

inline auto pretty_name(
    std::string& name,
    std::map<std::string, std::string>& greek_letter_to_latex) -> std::string& {
  for (const auto& item : greek_letter_to_latex) {
    if (name.find(item.first) != std::string::npos) {
      // Replace the Greek letter with its LaTeX representation
      // Lambda2 -> {\\lambda}2
      name =
          name.replace(name.find(item.first), item.first.length(), item.second);
    }
  }
  return name;
}

auto generate_markdown_table(
    const Settings& settings,
    const std::vector<ConstituentId>& modeled_constituents) -> std::string {
  std::map<std::string, std::string> greek_letter_to_latex = {
      {"Alpha", "{\\alpha}"}, {"Beta", "{\\beta}"}, {"Gamma", "{\\gamma}"},
      {"Delta", "{\\delta}"}, {"Psi", "{\\psi}"},   {"Phi", "{\\phi}"},
      {"Theta", "{\\theta}"}, {"Chi", "{\\chi}"},   {"Pi", "{\\pi}"},
      {"Mu", "{\\mu}"},       {"Nu", "{\\nu}"},     {"Lambda", "{\\lambda}"},
      {"Eps", "{\\epsilon}"}, {"Eta", "{\\eta}"},   {"Sigma", "{\\sigma}"},
      {"Ups", "{\\upsilon}"}, {"Rho", "{\\rho}"},   {"Tau", "{\\tau}"},
  };

  auto wt = wave_table_factory(settings.engine_type());
  wt->set_modeled_constituents(modeled_constituents);
  auto inference = inference_factory(*wt, settings.inference_type());

  auto inferred_constituents = inference->inferred_constituents();

  std::string table = "";

  table += "| Constituent | Speed (Deg/hr) | XDO | Modeled | Inferred |\n";
  table += "|-------------|----------------|-----|---------|----------|\n";
  for (const auto& item : *wt) {
    const auto& wave = item.value();
    const auto ident = wave->ident();
    auto name = std::string(wave->name());
    pretty_name(name, greek_letter_to_latex);
    const auto speed = wave->frequency<kDegreePerHour>();
    const auto xdo = wave->xdo_alphabetical();

    const auto is_modeled = wave->is_modeled() ? "Yes" : "No";
    const auto is_inferred =
        check_if_is_inferred(ident, wave->is_modeled(), inferred_constituents)
            ? "Yes"
            : "No";

    table += "| " + name + " | " + std::to_string(speed) + " | " + xdo + " | " +
             is_modeled + " | " + is_inferred + " |\n";
  }
  table += "\n";
  return table;
}

}  // namespace fes