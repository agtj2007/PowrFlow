function [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, Ybus, Yf, Yt, V, ref, pv, pq);
%PFSOLN  Updates bus, gen, branch data structures to match power flow soln.
%   [bus, gen, branch] = pfsoln(baseMVA, bus0, gen0, branch0, ...
%                                   Ybus, Yf, Yt, V, ref, pv, pq)

%   MATPOWER
%   $Id: pfsoln.m,v 1.6 2004/08/23 20:56:52 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% constants
j = sqrt(-1);
nl = size(branch0, 1);      %% number of lines

%% define named indices into bus, gen, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%% initialize return values
bus     = bus0;
gen     = gen0;
branch  = branch0;

%%----- update bus voltages -----
bus(:, VM) = abs(V);
bus(:, VA) = angle(V) * 180 / pi;

%%----- update Qg for all gens and Pg for swing bus -----
%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?
refgen = find(gbus == ref);             %% which is(are) the reference gen(s)?

%% compute total injected bus powers
%% This is slow in Matlab 5 ...
% Sg = V(gbus) .* conj(Ybus(gbus, :) * V);
%% ... so we do this instead ...
temp = Ybus.';
Sg = V(gbus) .* conj(temp(:, gbus).' * V);

%% update Qg for all generators
gen(:, QG) = zeros(size(gen, 1), 1);                %% zero out all Qg
gen(on, QG) = imag(Sg) * baseMVA + bus(gbus, QD);   %% inj Q + local Qd
%% ... at this point any buses with more than one generator will have
%% the total Q dispatch for the bus assigned to each generator. This
%% must be split between them. We do it equally.

if length(on) > 1
    %% build connection matrix, element i, j is 1 if gen on(j) at bus i is ON
    nb = size(bus, 1);
    ngon = size(on, 1);
    Cg = sparse(gbus, [1:ngon]', ones(ngon, 1), nb, ngon);

    %% divide Qg by number of generators at the bus to distribute equally
    gen(on, QG) = gen(on, QG) ./ (Cg' * sum(Cg')');
end

%% update Pg for swing bus
gen(on(refgen(1)), PG) = real(Sg(refgen(1))) * baseMVA + bus(ref, PD);  %% inj P + local Pd
if length(refgen) > 1       %% more than one generator at the ref bus
    %% subtract off what is generated by other gens at this bus
    gen(on(refgen(1)), PG) = gen(on(refgen(1)), PG) - sum(gen(on(refgen(2:length(refgen))), PG));
end

%%----- update/compute branch power flows -----
out = find(branch(:, BR_STATUS) == 0);      %% out-of-service branches
br = find(branch(:, BR_STATUS));            %% in-service branches
Sf = V(branch(br, F_BUS)) .* conj(Yf(br, :) * V) * baseMVA; %% complex power at "from" bus
St = V(branch(br, T_BUS)) .* conj(Yt(br, :) * V) * baseMVA; %% complex power injected at "to" bus
branch(br, [PF, QF, PT, QT]) = [real(Sf) imag(Sf) real(St) imag(St)];
branch(out, [PF, QF, PT, QT]) = zeros(length(out), 4);

return;
