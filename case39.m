function [baseMVA, bus, gen, branch, areas, gencost] = case39
%CASE39    Power flow data for 39 bus case.
%   Please see 'help caseformat' for details on the case file format.
%
%   Based on data from:
%
%       M. A. Pai, Energy Function Analysis for Power System Stability,
%       Kluwer Academic Publishers, Boston, 1989.

%   MATPOWER
%   $Id: case39.m,v 1.6 2004/09/21 01:48:33 ray Exp $

%%-----  Power Flow Data  -----%%
%% system MVA base
baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
bus = [
	1	1	0	0	0	0	1	1	0	0	1	1.06	0.94;
	2	1	0	0	0	0	1	1	0	0	1	1.06	0.94;
	3	1	322	2.4	0	0	1	1.0341	-9.73	0	1	1.06	0.94;
	4	1	500	184	0	0	1	1.0116	-10.53	0	1	1.06	0.94;
	5	1	0	0	0	0	1	1.0165	-9.38	0	1	1.06	0.94;
	6	1	0	0	0	0	1	1.0172	-8.68	0	1	1.06	0.94;
	7	1	233.8	84	0	0	1	1.0067	-10.84	0	1	1.06	0.94;
	8	1	522	176.6	0	0	1	1.0057	-11.34	0	1	1.06	0.94;
	9	1	0	0	0	0	1	1.0322	-11.15	0	1	1.06	0.94;
	10	1	0	0	0	0	1	1.0235	-6.31	0	1	1.06	0.94;
	11	1	0	0	0	0	1	1.0201	-7.12	0	1	1.06	0.94;
	12	1	8.5	88	0	0	1	1.0072	-7.14	0	1	1.06	0.94;
	13	1	0	0	0	0	1	1.0207	-7.02	0	1	1.06	0.94;
	14	1	0	0	0	0	1	1.0181	-8.66	0	1	1.06	0.94;
	15	1	320	153	0	0	1	1.0194	-9.06	0	1	1.06	0.94;
	16	1	329.4	32.3	0	0	1	1.0346	-7.66	0	1	1.06	0.94;
	17	1	0	0	0	0	1	1.0365	-8.65	0	1	1.06	0.94;
	18	1	158	30	0	0	1	1.0343	-9.49	0	1	1.06	0.94;
	19	1	0	0	0	0	1	1.0509	-3.04	0	1	1.06	0.94;
	20	1	680	103	0	0	1	0.9914	-4.45	0	1	1.06	0.94;
	21	1	274	115	0	0	1	1.0337	-5.26	0	1	1.06	0.94;
	22	1	0	0	0	0	1	1.0509	-0.82	0	1	1.06	0.94;
	23	1	247.5	84.6	0	0	1	1.0459	-1.02	0	1	1.06	0.94;
	24	1	308.6	-92.2	0	0	1	1.0399	-7.54	0	1	1.06	0.94;
	25	1	224	47.2	0	0	1	1.0587	-5.51	0	1	1.06	0.94;
	26	1	139	17	0	0	1	1.0536	-6.77	0	1	1.06	0.94;
	27	1	281	75.5	0	0	1	1.0399	-8.78	0	1	1.06	0.94;
	28	1	206	27.6	0	0	1	1.0509	-3.27	0	1	1.06	0.94;
	29	1	283.5	26.9	0	0	1	1.0505	-0.51	0	1	1.06	0.94;
	30	2	0	0	0	0	1	1.0475	0	0	1	1.06	0.94;
	31	3	9.2	4.6	0	0	1	0.982	0	0	1	1.06	0.94;
	32	2	0	0	0	0	1	0.9831	1.63	0	1	1.06	0.94;
	33	2	0	0	0	0	1	0.9972	2.18	0	1	1.06	0.94;
	34	2	0	0	0	0	1	1.0123	0.74	0	1	1.06	0.94;
	35	2	0	0	0	0	1	1.0493	4.14	0	1	1.06	0.94;
	36	2	0	0	0	0	1	1.0635	6.83	0	1	1.06	0.94;
	37	2	0	0	0	0	1	1.0278	1.27	0	1	1.06	0.94;
	38	2	0	0	0	0	1	1.0265	6.55	0	1	1.06	0.94;
	39	2	1104	250	0	0	1	1.03	-10.96	0	1	1.06	0.94;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
gen = [
	30	250	103.3	9999	-9999	1.0475	100	1	350	0;
	31	572.9	170.3	9999	-9999	0.982	100	1	1145.55	0;
	32	650	175.9	9999	-9999	0.9831	100	1	750	0;
	33	632	103.3	9999	-9999	0.9972	100	1	732	0;
	34	508	164.4	9999	-9999	1.0123	100	1	608	0;
	35	650	204.8	9999	-9999	1.0493	100	1	750	0;
	36	560	96.9	9999	-9999	1.0635	100	1	660	0;
	37	540	-4.4	9999	-9999	1.0278	100	1	640	0;
	38	830	19.4	9999	-9999	1.0265	100	1	930	0;
	39	1000	68.5	9999	-9999	1.03	100	1	1100	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status
branch = [
	1	2	0.0035	0.0411	0.6987	9900	0	0	0	0	1;
	1	39	0.001	0.025	0.75	9900	0	0	0	0	1;
	2	3	0.0013	0.0151	0.2572	9900	0	0	0	0	1;
	2	25	0.007	0.0086	0.146	9900	0	0	0	0	1;
	3	4	0.0013	0.0213	0.2214	9900	0	0	0	0	1;
	3	18	0.0011	0.0133	0.2138	9900	0	0	0	0	1;
	4	5	0.0008	0.0128	0.1342	9900	0	0	0	0	1;
	4	14	0.0008	0.0129	0.1382	9900	0	0	0	0	1;
	5	6	0.0002	0.0026	0.0434	9900	0	0	0	0	1;
	5	8	0.0008	0.0112	0.1476	9900	0	0	0	0	1;
	6	7	0.0006	0.0092	0.113	9900	0	0	0	0	1;
	6	11	0.0007	0.0082	0.1389	9900	0	0	0	0	1;
	7	8	0.0004	0.0046	0.078	9900	0	0	0	0	1;
	8	9	0.0023	0.0363	0.3804	9900	0	0	0	0	1;
	9	39	0.001	0.025	1.2	9900	0	0	0	0	1;
	10	11	0.0004	0.0043	0.0729	9900	0	0	0	0	1;
	10	13	0.0004	0.0043	0.0729	9900	0	0	0	0	1;
	13	14	0.0009	0.0101	0.1723	9900	0	0	0	0	1;
	14	15	0.0018	0.0217	0.366	9900	0	0	0	0	1;
	15	16	0.0009	0.0094	0.171	9900	0	0	0	0	1;
	16	17	0.0007	0.0089	0.1342	9900	0	0	0	0	1;
	16	19	0.0016	0.0195	0.304	9900	0	0	0	0	1;
	16	21	0.0008	0.0135	0.2548	9900	0	0	0	0	1;
	16	24	0.0003	0.0059	0.068	9900	0	0	0	0	1;
	17	18	0.0007	0.0082	0.1319	9900	0	0	0	0	1;
	17	27	0.0013	0.0173	0.3216	9900	0	0	0	0	1;
	21	22	0.0008	0.014	0.2565	9900	0	0	0	0	1;
	22	23	0.0006	0.0096	0.1846	9900	0	0	0	0	1;
	23	24	0.0022	0.035	0.361	9900	0	0	0	0	1;
	25	26	0.0032	0.0323	0.513	9900	0	0	0	0	1;
	26	27	0.0014	0.0147	0.2396	9900	0	0	0	0	1;
	26	28	0.0043	0.0474	0.7802	9900	0	0	0	0	1;
	26	29	0.0057	0.0625	1.029	9900	0	0	0	0	1;
	28	29	0.0014	0.0151	0.249	9900	0	0	0	0	1;
	12	11	0.0016	0.0435	0	9900	0	0	1.006	0	1;
	12	13	0.0016	0.0435	0	9900	0	0	1.006	0	1;
	6	31	0	0.025	0	9900	0	0	1.07	0	1;
	10	32	0	0.02	0	9900	0	0	1.07	0	1;
	19	33	0.0007	0.0142	0	9900	0	0	1.07	0	1;
	20	34	0.0009	0.018	0	9900	0	0	1.009	0	1;
	22	35	0	0.0143	0	9900	0	0	1.025	0	1;
	23	36	0.0005	0.0272	0	9900	0	0	1	0	1;
	25	37	0.0006	0.0232	0	9900	0	0	1.025	0	1;
	2	30	0	0.0181	0	9900	0	0	1.025	0	1;
	29	38	0.0008	0.0156	0	9900	0	0	1.025	0	1;
	19	20	0.0007	0.0138	0	9900	0	0	1.06	0	1;
];

%%-----  OPF Data  -----%%
%% area data
areas = [
	1	1;
];

%% generator cost data
%	1	startup	shutdown	n	x0	y0	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
gencost = [
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.01	0.3	0.2;
	2	0	0	3	0.006	0.3	0.2;
	2	0	0	3	0.006	0.3	0.2;
];

return;
