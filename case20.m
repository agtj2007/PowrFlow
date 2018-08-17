function [baseMVA, bus, gen, branch, areas, gencost] = case20
%CASE30    Power flow data for 30 bus, 6 generator case.
%   Please see 'help caseformat' for details on the case file format.
%
%   Based on data from ...
%     Alsac, O. & Stott, B., "Optimal Load Flow with Steady State Security",
%     IEEE Transactions on Power Apparatus and Systems, Vol. PAS 93, No. 3,
%     1974, pp. 745-751.
%   ... with branch parameters rounded to nearest 0.01, shunt values divided
%   by 100 and shunt on bus 10 moved to bus 5, load at bus 5 zeroed out.
%   Generator locations, costs and limits and bus areas were taken from ...
%     Ferrero, R.W., Shahidehpour, S.M., Ramesh, V.C., "Transaction analysis
%     in deregulated power systems using game theory", IEEE Transactions on
%     Power Systems, Vol. 12, No. 3, Aug 1997, pp. 1340-1347.
%   Generator Q limits were derived from Alsac & Stott, using their Pmax
%   capacities. V limits and line |S| limits taken from Alsac & Stott.

%   MATPOWER
%   $Id: case30.m,v 1.7 2004/09/21 02:41:49 ray Exp $

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 10;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	3	0	0	0	0	1	1.02	0	10	0	1.05	0.95;
	2	2	0.0497	0.0261	0	0	1	1	0	10	1	1.1	0.95;
	3	2	0.0034	0.0065	0	0	1	1	0	10	1	1.05	0.95;
	4	2	0.0487	0.0002	0	0	1	1	0	10	1	1.05	0.95;
	5	2	0.0461	0.0163	0	0.19	1	1	0	10	1	1.05	0.95;
	6	2	0.0479	0.0180	0	0	1	1	0	10	1	1.05	0.95;
	7	2	0.0282	0.0497	0	0	1	1	0	10	1	1.05	0.95;
	8	1	0.0475	0.0470	0	0	1	1	0	10	1	1.05	0.95;
	9	2	0.0159	0.0217	0	0	1	1	0	10	1	1.05	0.95;
	10	2	0.0133	0.0468	0	0	3	1	0	10	1	1.05	0.95;
	11	1	0.0031	0.0248	0	0	1	1	0	10	1	1.05	0.95;
	12	2	0.0166	0.0486	0	0	2	1	0	10	1	1.05	0.95;
	13	1	0.0146	0.0479	0	0	2	1	0	10	1	1.1	0.95;
	14	1	0.0201	0.0360	0	0	2	1	0	10	1	1.05	0.95;
	15	2	0.0017	0.0324	0	0	2	1	0	10	1	1.05	0.95;
	16	2	0.0144	0.0431	0	0	2	1	0	10	1	1.05	0.95;
	17	2	0.0110	0.0015	0	0	2	1	0	10	1	1.05	0.95;
	18	2	0.0283	0.0183	0	0	2	1	0	10	1	1.05	0.95;
	19	2	0.0435	0.0047	0	0	2	1	0	10	1	1.05	0.95;
	20	2	0.0127	0.0167	0	0	2	1	0	10	1	1.05	0.95;
	21	2	0.0296	0.0004	0	0	3	1	0	10	1	1.05	0.95;
	];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
gen = [
	1	0.0200	0.0018	0.1	-0.1	1	10	1	80	0;
    5	0.0025	0.0018	0.1	-0.1	1	10	1	80	0;
	8	0.0034	0.0025	0.1	-0.1	1	100	1	80	0;
	11	0.0424	0.0309	0.1	-0.1    1	100	1	50	0;
	13	0.0677  0.0505	0.1	-0.1	1	100	1	55	0;
	14	0.0816	0.0598	0.1 -0.1	1	100	1	30	0;	
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
branch = [
	1	2	0.0002	0.0005	0.03	130	130	130	0	0	1;
	2	3	0.0082	0.0246	0.02	130	130	130	0	0	1;
	3	4	0.0006	0.0017	0.02	65	65	65	0	0	1;
	4	5	0.0009	0.0027	0	130	130	130	0	0	1;
	5	6	0.0070	0.0209	0.02	130	130	130	0	0	1;
	6	7	0.0015	0.0044	0.02	65	65	65	0	0	1;
	7	8	0.0006	0.0017	0	90	90	90	0	0	1;
    8	9	0.0006	0.0019	0	90	90	90	0	0	1;
	9	10	0.0072	0.0217	0.01	70	70	70	0	0	1;
	10	11	0.0003	0.0010	0.01	130	130	130	0	0	1;
	11	12	0.0070	0.0209	0	32	32	32	0	0	1;
	12	13	0.0022	0.0066	0	65	65	65	0	0	1;
	13	14	0.0063	0.0188	0	32	32	32	0	0	1;
	14	15	0.0051	0.0154	0	65	65	65	0	0	1;
	15	16	0.0058	0.0173	0	65	65	65	0	0	1;
	16	17	0.0029	0.0086	0	65	65	65	0	0	1;
	17	18	0.0089	0.0268	0	65	65	65	0	0	1;
	18	19	0.0009	0.0028	0	32	32	32	0	0	1;
	19	20	0.0029	0.0088	0	32	32	32	0	0	1;
	20	21	0.0082	0.0247	0	32	32	32	0	0	1;	
];

%%-----  OPF Data  -----%%
%% area data
areas = [
	1	8;
	2	23;
	3	26;
];

%% generator cost data
%	1	startup	shutdown	n	x0	y0	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.02	2	0;
	2	0	0	3	0.0175	1.75	0;
	2	0	0	3	0.0625	1	0;
	2	0	0	3	0.00834	3.25	0;
	2	0	0	3	0.025	3	0;
	2	0	0	3	0.025	3	0;
];

return;
