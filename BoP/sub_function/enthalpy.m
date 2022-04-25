function hf = enthalpy(id, T)
% The function enthalpy calculates enthalpy of formation of relevant
% chemicals in SOFC system. Taken from https://webbook.nist.gov/chemistry/
% Arguments:
%   _id: chemical molecule short hand. Special case: H2Ol is liquid water
%   _T: temperature in Kelvin
% Return:
%   _hf: enthalpy of formation in J/mol. No return if given 'id' was
%   invalid
switch(id)
    case 'H2'
        if T<1000 %298-1000. K
            A=33.066178;B=-11.363417;C=11.432816;D=-2.772874;E=-0.158558;F=-9.980797;
        else %1000. - 2500. K
            A=18.563083;B=12.257357;C=-2.859786;D=0.268238;E=1.97799;F=-1.147438;
        end
    case 'H2O' %500. - 1700. K, gas phase
        A=30.092;B=6.832514;C=6.793435;D=-2.53448;E=0.082139;F=-250.881;
    case 'H2Ol' %298. - 500. K, condensed phase
        A=-203.6060;B=1523.290;C=-3196.413;D=2474.455;E=3.855326;F=-256.5478;
    case 'O2'
        if T<700 %100. - 700.K
            A=31.32234;B=-20.23531;C=57.86644;D=-36.50624;E=-0.007374;F=-8.903471;
        else % 700. - 2000.	K
            A=30.03235;B=8.772972;C=-3.988133;D=0.788313;E=-0.741599;F=-11.32468;
        end
    case 'N2'
        if T<500 %100. - 500. K
            A=28.98641;B=1.853978;C=-9.647459;D=16.63537;E=0.000117;F=-8.671914;
        else %500. - 2000.K
            A=19.50583;B=19.88705;C=-8.598535;D=1.369784;E=0.527601;F=-4.935202;
        end
    case 'CH4' %298. - 1300. K
        A=-0.703029;B=108.4773;C=-42.52157;D=5.862788;E=0.678565;F=-76.84376;
    case 'CO' %298. - 1300.	
        A=25.56759;B=6.09613;C=4.054656;D=-2.671301;E=0.131021;F=-118.0089;
    case 'CO2' %298. - 1200. K
        A=24.99735;B=55.18696;C=-33.69137;D=7.948387;E=-0.136638;F=-403.6075;
    case 'SMR' % CH4 + H2O == CO + 3H2
        hf = enthalpy("CO",T) + 3*enthalpy("H2",T) - enthalpy("H2O",T) - enthalpy("CH4",T);
        return;
    case 'WGS' % CO + H2O == CO2 + H2
        hf = enthalpy("CO2",T) + enthalpy("H2",T) - enthalpy("CO",T) - enthalpy("H2O",T);
        return;
    case 'EL'
        hf = enthalpy("H2O",T) - enthalpy("H2",T) - 1/2*enthalpy("O2",T);
        return;
    case 'evap'
        hf = enthalpy("H2O",T) - enthalpy("H2Ol",T);
        return;
    otherwise
        fprintf('invalid ID \n');
        return;
end

%find enthalpy at T [K]:
t = T/1000;
hf = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F;
hf = hf*1000;
end
