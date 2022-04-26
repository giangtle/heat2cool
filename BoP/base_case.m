function [CH4, H2O, air] = base_case(T_env, H_env)
% This function returns the base case feed streams to the fuel cell
%   _ natural gas flow has a LHV of 2500 watts.
%   _ water stream gives a steam to gas ratio of 2:1.
%   _ the air flow gives an air to gas ratio of 35:1.

% Natural gas stream to FC:
CH4 = struct(); CH4.phase = "gas"; CH4.n = SLM_to_mol_sec(3.954); CH4.yCH4 = 1; CH4.T = T_env;
% Water stream to FC based on 2:1 steam to gas ratio:
H2O = struct(); H2O.phase = "liq"; H2O.yH2Ol = 1; H2O.T = T_env; H2O.n = CH4.n*2;
% Air stream to FC:
air = struct(); air.phase = "gas"; air.n = CH4.n*35; air.T = T_env;
air.yH2O = H_env*cal_yH2Osat(T_env); air.yO2 = 0.21*(1-air.yH2O); air.yN2 = 0.79*(1-air.yH2O);
end
