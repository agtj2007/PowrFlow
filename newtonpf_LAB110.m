function [V, converged, i] = newtonpf_LB110(Ybus, Sbus, V0, ref, pv, pq, nb,mpopt)
; 
%NEWTONPF  Solves the power flow using a full Newton's method.
%   [V, converged, i] = newtonpf(Ybus, Sbus, V0, ref, pv, pq, mpopt)
%   solves for bus voltages given the full system admittance matrix (for
%   all buses), the complex bus power injection vector (for all buses),
%   the initial vector of complex bus voltages, and column vectors with
%   the lists of bus indices for the swing bus, PV buses, and PQ buses,
%   respectively. The bus voltage vector contains the set point for
%   generator (including ref bus) buses, and the reference angle of the
%   swing bus, as well as an initial guess for remaining magnitudes and
%   angles. mpopt is a MATPOWER options vector which can be used to 
%   set the termination tolerance, maximum number of iterations, and 
%   output options (see 'help mpoption' for details). Uses default
%   options if this parameter is not given. Returns the final complex
%   voltages, a flag which indicates whether it converged or not, and
%   the number of iterations performed.

%   MATPOWER
%   $Id: newtonpf.m,v 1.5 2004/08/23 20:56:46 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% default arguments
if nargin < 7
    mpopt = mpoption;
end

%% options
tol     = mpopt(2);
max_it  = mpopt(3);
verbose = mpopt(31);

%% initialize
j = sqrt(-1);
converged = 0;
counter = 0;
V = V0;
Va = angle(V);
Vm = abs(V);
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

%% evaluate F(x0)
P=zeros(nb,1);
Q=zeros(nb,1);
for i=1:nb
    for k=1:nb 
        if G(i,k)~=0 ||B(i,k)~=0   %zero elements do not take part in calculation. A good mesure to speed up calculation. 
                                     % Almost as fast as calculation in matpower.
            P(i)=P(i)+Vm(i)*Vm(k)*(G(i,k)*cos(Va(i)-Va(k))+B(i,k)*sin(Va(i)-Va(k)));      
            Q(i)=Q(i)+Vm(i)*Vm(k)*(G(i,k)*sin(Va(i)-Va(k))-B(i,k)*cos(Va(i)-Va(k)));
        end
    end      
end
for i=1:nb  
        dP(i)=P(i)-Ps(i);   
        dQ(i)=Q(i)-Qs(i);  
end
F = [  dP([pv;pq])';
        dQ([pq])'];
%% check tolerance
normF = norm(F, inf);
if verbose > 1
    fprintf('\n it    max P & Q mismatch (p.u.)');
    fprintf('\n----  ---------------------------');
    fprintf('\n%3d        %10.3e', i, normF);
end
if normF < tol
    converged = 1;
    if verbose > 1
        fprintf('\nConverged!\n');
    end
end

%% do Newton iterations
while (~converged & counter < max_it)
    %% update iteration counter
    counter = counter + 1;    
    J=JaccobiMatrix_LAB110(Va,Vm,nb,npq,npv,G,B,P,Q,pv,pq,ref); 
    
    %% compute update step
    dx = -(J \ F);    

    %% update voltage
    if npv
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    if npq
        Va(pq) = Va(pq) + dx(j3:j4);
        Vm(pq) = Vm(pq) + dx(j5:j6);
    end
    V = Vm .* exp(j * Va);
    Vm = abs(V);            %% update Vm and Va again in case
    Va = angle(V);          %% we wrapped around with a negative Vm

    %% evalute F(x)
    P=zeros(nb,1);
    Q=zeros(nb,1);
    for i=1:nb
        for k=1:nb        
            if G(i,k)~=0 ||B(i,k)~=0  % zero elements do not take part in calculation. A good mesure to speed up calculation. 
                                     % Almost as fast as calculation in matpower.
                  P(i)=P(i)+Vm(i)*Vm(k)*(G(i,k)*cos(Va(i)-Va(k))+B(i,k)*sin(Va(i)-Va(k)));      
                  Q(i)=Q(i)+Vm(i)*Vm(k)*(G(i,k)*sin(Va(i)-Va(k))-B(i,k)*cos(Va(i)-Va(k)));
            end
        end
    end
    for i=1:nb  
        dP(i)=P(i)-Ps(i);   
        dQ(i)=Q(i)-Qs(i);  
    end
    F = [  dP([pv;pq])';
        dQ([pq])'];

    %% check for convergence
    normF = norm(F, inf);
    if verbose > 1
        fprintf('\n%3d        %10.3e', counter, normF);
    end
    if normF < tol
        converged = 1;
        if verbose
            fprintf('\nNewton''s method power flow converged in %d iterations.\n', counter);
        end
    end
end

if verbose
    if ~converged
        fprintf('\nNewton''s method power did not converge in %d iterations.\n', counter);
    end
end
