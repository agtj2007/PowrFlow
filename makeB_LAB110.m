function [Bp, Bpp] = makeB_LAB110(baseMVA, bus, branch, alg)
%MAKEB   Builds the FDPF matrices, B prime and B double prime.
%   [Bp, Bpp] = makeB(baseMVA, bus, branch, alg) returns the two
%   matrices B prime and B double prime used in the fast decoupled power
%   flow. Does appropriate conversions to p.u.

%   MATPOWER
%   $Id: makeB.m,v 1.4 2004/08/23 20:56:39 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%   �滻���Bp��Bpp����� 2016-10-21

%% constants
nb = size(bus, 1);          %% number of buses
nl = size(branch, 1);       %% number of lines

%% define named indices into bus, branch matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;

%%-----  form Bp (B prime)  -----
temp_branch = branch;                       %% modify a copy of branch
temp_bus = bus;                           %% modify a copy of bus
for i=1:nb
    temp_bus(i,BS)=0;                      %% zero out shunts at buses
end
for k=1:nl
    temp_branch(k,BR_B)=0;                 %% zero out line charging shunts 
    temp_branch(k,TAP)=1;                  %% cancel out taps 
    if alg==2
        temp_branch(k,BR_R)=0;             %% zero out line resistance
    end
end
Bp = -imag( makeYbus_LAB110(baseMVA, temp_bus, temp_branch));
%%-----  form Bpp (B double prime)  -----
if nargout == 2
    temp_branch = branch;                       %% modify a copy of branch     
    for k=1:nl
        temp_branch(i,SHIFT)=0;                %% zero out phase shifters
        if alg==3                              %% if BX method
            temp_branch(i,BR_R)=0;             %% zero out line resistance
        end
    end
    Bpp = -imag( makeYbus_LAB110(baseMVA, bus, temp_branch) );
end
return;
