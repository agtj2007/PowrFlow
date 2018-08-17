function [V, converged, counter] = fdpf_LAB110(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%FDPF  Solves the power flow using a fast decoupled method.
%   [V, converged, i] = fdpf(Ybus, Sbus, V0, Bp, Bpp, ref, pv, pq, mpopt)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, the FDPF matrices B prime
%   and B double prime, and column vectors with the lists of bus indices
%   for the swing bus, PV buses, and PQ buses, respectively. The bus voltage
%   vector contains the set point for generator (including ref bus)
%   buses, and the reference angle of the swing bus, as well as an initial
%   guess for remaining magnitudes and angles. mpopt is a MATPOWER options
%   vector which can be used to set the termination tolerance, maximum
%   number of iterations, and  output options (see 'help mpoption'
%   for details). Uses default options if this parameter is not given.
%   Returns the final complex voltages, a flag which indicates whether it
%   converged or not, and the number of iterations performed.

%   MATPOWER
%   $Id: fdpf.m,v 1.5 2004/08/23 20:56:18 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%    
%% default arguments
if nargin < 7
    mpopt = mpoption;
end
%% options
tol     = mpopt(2);
max_it  = mpopt(4);
verbose = mpopt(31);
%% initialize
j = sqrt(-1);
converged = 0;
counter = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
nb=size(pv,1)+size(pq,1)+1;
G=real(Ybus);
B=imag(Ybus);
Ps=real(Sbus);
Qs=imag(Sbus);

%% set up indexing for updating V
npv = length(pv);
npq = length(pq);
j1 = 1;         j2 = npv;           %% j1:j2 - V angle of pv buses
j3 = j2 + 1;    j4 = j2 + npq;      %% j3:j4 - V angle of pq buses
j5 = j4 + 1;    j6 = j4 + npq;      %% j5:j6 - V mag of pq buses
%% evaluate initial mismatch
P=zeros(nb,1);
Q=zeros(nb,1);
dP=zeros(nb,1);
dQ=zeros(nb,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %node method, easy to understand, but slow
for i=1:nb                     
    for k=1:nb       
        if G(i,k)~=0 || B(i,k)~=0; % zero elements do not take part in calculation. A good mesure to speed up calculation. 
                                     % Almost as fast as calculation in matpower.
            P(i)=(P(i)+Vm(i)*Vm(k)*(G(i,k)*cos(Va(i)-Va(k))+B(i,k)*sin(Va(i)-Va(k))));      
            Q(i)=(Q(i)+Vm(i)*Vm(k)*(G(i,k)*sin(Va(i)-Va(k))-B(i,k)*cos(Va(i)-Va(k))));
        end
    end
    dP(i)=(P(i)-Ps(i))/Vm(i);   
    dQ(i)=(Q(i)-Qs(i))/Vm(i);  
end
PP=dP([pv;pq]);
QQ=dQ(pq);
%node method, easy to understand, but slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vector method, very fast in speed, but need more time to understand.
% mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;   
% P = real(mis([pv; pq]));
% Q = imag(mis(pq));
%vector method, very fast in speed, but need more time to understand.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check tolerance
normP = norm(PP, inf);
normQ = norm(QQ, inf);
if verbose > 1
    fprintf('\niteration     max mismatch (p.u.)  ');
    fprintf('\ntype   #        P            Q     ');
    fprintf('\n---- ----  -----------  -----------');
    fprintf('\n  -  %3d   %10.3e   %10.3e', i, normP, normQ);
end
if normP < tol && normQ < tol
    converged = 1;
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end
%% reduce B matrices
temp = Bp(:, [pv; pq])';
Bp = temp(:, [pv; pq])';
temp = Bpp(:, pq)';
Bpp = temp(:, pq)';
%% factor B matrices
[Lp, Up, Pp] = lu(Bp);
[Lpp, Upp, Ppp] = lu(Bpp);
%% do P and Q iterations
while (~converged && counter < max_it)
    %% update iteration counter
    counter = counter + 1;
    %%-----  do P iteration, update Va  -----
    dVa = -( Up \  (Lp \ (Pp * PP)));
    %% update voltage
    Va([pv; pq]) = Va([pv; pq]) + dVa;
    V = Vm .* exp(j * Va);
    Vm = abs(V);
  %% evalute mismatch
    P=zeros(nb,1);
    Q=zeros(nb,1);
    dP=zeros(nb,1);
    dQ=zeros(nb,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %node method
    for i=1:nb
        for k=1:nb
            if G(i,k)~=0 || B(i,k)~=0; % zero elements do not take part in calculation. A good mesure to speed up calculation.
             % Almost as fast as calculation in matpower.
                 P(i)=(P(i)+Vm(i)*Vm(k)*(G(i,k)*cos(Va(i)-Va(k))+B(i,k)*sin(Va(i)-Va(k))));
                 Q(i)=(Q(i)+Vm(i)*Vm(k)*(G(i,k)*sin(Va(i)-Va(k))-B(i,k)*cos(Va(i)-Va(k))));
            end
        end
        dP(i)=(P(i)-Ps(i))/Vm(i);
        dQ(i)=(Q(i)-Qs(i))/Vm(i);
    end     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %vector method, very fast in speed, but need more time to understand.
%     mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;   
%     P = real(mis([pv; pq]));
%     Q = imag(mis(pq));
    %vector method, very fast in speed, but need more time to understand.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PP=dP([pv;pq]);
    QQ=dQ(pq);
%% check tolerance
    normP = norm(PP, inf);
    normQ = norm(QQ, inf);
    if verbose > 1
        fprintf('\n  P  %3d   %10.3e   %10.3e', i, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n', counter, counter-1);
        end
        break;
    end
    %%-----  do Q iteration, update Vm  -----
    dVm = -( Upp \ (Lpp \ (Ppp * QQ)) );
    %% update voltage
    Vm(pq) = Vm(pq) + dVm;
    V = Vm .* exp(j * Va);
    Va = angle(V);
    %% evalute mismatch
    P=zeros(nb,1);
    Q=zeros(nb,1);
    dP=zeros(nb,1);
    dQ=zeros(nb,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %node method
    for i=1:nb                     
        for k=1:nb      
            if G(i,k)~=0 || B(i,k)~=0;   % zero elements of G and B matrix do not take part in calculation. A good mesure to speed up calculation. 
                                     % Almost as fast as calculation in matpower.
                 P(i)=(P(i)+Vm(i)*Vm(k)*(G(i,k)*cos(Va(i)-Va(k))+B(i,k)*sin(Va(i)-Va(k))));      
                 Q(i)=(Q(i)+Vm(i)*Vm(k)*(G(i,k)*sin(Va(i)-Va(k))-B(i,k)*cos(Va(i)-Va(k))));
            end
        end
        dP(i)=(P(i)-Ps(i))/Vm(i);   
        dQ(i)=(Q(i)-Qs(i))/Vm(i);  
        end
%node method, easy to understand, but slow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %vector method, very fast in speed, but need more time to understand.
%       mis = (V .* conj(Ybus * V) - Sbus) ./ Vm;   
%       P = real(mis([pv; pq]));
%       Q = imag(mis(pq));
     %vector method, very fast in speed, but need more time to understand.
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PP=dP([pv;pq]);
    QQ=dQ(pq);    
    normP = norm(PP, inf);
    normQ = norm(QQ, inf);          
    if verbose > 1
        fprintf('\n  Q  %3d   %10.3e   %10.3e', counter, normP, normQ);
    end
    if normP < tol && normQ < tol
        converged = 1;
        if verbose
            fprintf('\nFast-decoupled power flow converged in %d P-iterations and %d Q-iterations.\n',counter,counter);
        end
        break;
    end
end
if verbose
    if ~converged
        fprintf('\nFast-decoupled power flow did not converge in %d iterations.\n', counter);
    end
end
