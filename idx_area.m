function [AREA_I, PRICE_REF_BUS] = idx_area
%IDX_AREA   Defines variables for column indices to areas.
%   [AREA_I, PRICE_REF_BUS] = idx_area

%   MATPOWER
%   $Id: idx_area.m,v 1.5 2004/08/23 20:56:26 ray Exp $
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2004 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% define the indices
AREA_I          = 1;    %% area number
PRICE_REF_BUS   = 2;    %% price reference bus for this area

return;
