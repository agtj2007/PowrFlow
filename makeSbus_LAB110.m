function Sbus = makeSbus(baseMVA, bus, gen)
%MAKESBUS   Builds the vector of complex bus power injections.
%   Sbus = makeSbus(baseMVA, bus, gen) returns the vector of complex bus
%   power injections, that is, generation minus load. Power is expressed
%   in per unit.

%   MATPOWER
%   $Id: makeSbus.m,v 1.4 2004/08/23 20:56:41 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%   �޸�˵�����ð��ڵ����ע�빦�ʵķ������ԭ���������㷽����2016-10-21 
%% constants
j = sqrt(-1);

%% define named indices into bus, gen matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
    GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;

%% generator info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

nbus=size(bus,1);
ng=size(gen,1);
Pbs=zeros(nbus,1);
Qbs=zeros(nbus,1);
Pgs=zeros(nbus,1);
Qgs=zeros(nbus,1);
Ps=zeros(nbus,1);
Qs=zeros(nbus,1);
S=zeros(nbus,1);

 for iGen=1:ng               %���ҿ���ķ��������ĸ�߱��
     if (gen(iGen,GEN_STATUS)~=0)
         Busg(iGen)=gen(iGen,GEN_BUS);
     end
 end    
             
for iGbus=Busg                      %����������ӳ�䵽����нڵ���Ŀ�Ĺ��������ϣ��������ע�빦�ʼ���
    for i=1:ng
        if gen(i,GEN_BUS)==iGbus
            Pgs(iGbus)=gen(i,PG);
            Qgs(iGbus)=gen(i,QG);
        end
    end  
end

for iBus=1:nbus                     %ȡ���ڵ�ĸ���
    Pbs(iBus)=bus(iBus,PD);
    Qbs(iBus)=bus(iBus,QD);    
end

    Sbus=(Pgs-Pbs+j*(Qgs-Qbs))/baseMVA;%% power injected by generators
%% plus power injected by loads converted to p.u.   
return;
