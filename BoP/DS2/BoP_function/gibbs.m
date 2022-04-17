function gf = gibbs(id, T)
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
            A=33.066178;B=-11.363417;C=11.432816;D=-2.772874;E=-0.158558;F=-9.980797;H=0;G=172.707974;
        else %1000. - 2500. K
            A=18.563083;B=12.257357;C=-2.859786;D=0.268238;E=1.97799;F=-1.147438;H=0;G=156.288133;
        end
    case 'H2O' %500. - 1700. K, gas phase
        A=30.092;B=6.832514;C=6.793435;D=-2.53448;E=0.082139;F=-250.881;H=0;G=223.3967;
    case 'H2Ol' %298. - 500. K, condensed phase
        A=-203.6060;B=1523.290;C=-3196.413;D=2474.455;E=3.855326;F=-256.5478;H=0;G=-488.7163;
    case 'O2'
        if T<700 %100. - 700.K
            A=31.32234;B=-20.23531;C=57.86644;D=-36.50624;E=-0.007374;F=-8.903471;H=0;G=246.7945;
        else % 700. - 2000.	K
            A=30.03235;B=8.772972;C=-3.988133;D=0.788313;E=-0.741599;F=-11.32468;H=0;G=236.1663;
        end
    case 'N2'
        if T<500 %100. - 500. K
            A=28.98641;B=1.853978;C=-9.647459;D=16.63537;E=0.000117;F=-8.671914;H=0;G=226.4168;
        else %500. - 2000.K
            A=19.50583;B=19.88705;C=-8.598535;D=1.369784;E=0.527601;F=-4.935202;H=0;G=212.3900;
        end
    case 'CH4' %298. - 1300. K
        A=-0.703029;B=108.4773;C=-42.52157;D=5.862788;E=0.678565;F=-76.84376;H=0;G=158.7163;
    case 'CO' %298. - 1300.	
        A=25.56759;B=6.09613;C=4.054656;D=-2.671301;E=0.131021;F=-118.0089;H=0;G=227.3665;
    case 'CO2' %298. - 1200. K
        A=24.99735;B=55.18696;C=-33.69137;D=7.948387;E=-0.136638;F=-403.6075;H=0;G=228.2431;
    case 'SMR' % CH4 + H2O == CO + 3H2
        gf = gibbs("CO",T) + 3*gibbs("H2",T) - gibbs("H2O",T) - gibbs("CH4",T);
        return;
    case 'WGS' % CO + H2O == CO2 + H2
        gf = gibbs("CO2",T) + gibbs("H2",T) - gibbs("CO",T) - gibbs("H2O",T);
        return;
    case 'EL'
        gf = gibbs("H2O",T) - gibbs("H2",T) - 1/2*gibbs("O2",T);
        return;
    case 'evap'
        gf = gibbs("H2O",T) - gibbs("H2Ol",T);
        return;
    otherwise
        fprintf('invalid ID \n');
        return;
end

%find enthalpy at T [K]:
t = T/1000;
hf = A*t + B*t.^2/2 + C*t.^3/3 + D*t.^4/4 - E./t + F - H;
hf = hf*1000;
s = A*log(t) + B*t + C*t.^2/2 + D*t.^3/3 - E./(2*t.^2) + G;
gf = hf - T.*s;
end