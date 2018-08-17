       
function [Ybus,Yt,Yf] = makeYbus_LB110(baseMVA,bus,branch)
[nb,mb]=size(bus);
[nl,ml]=size(branch);
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
%      
Ybr=zeros(nb,nb);      
Ygb=zeros(nb,nb);
Ybus=zeros(nb,nb);% 
Ytt=zeros(nb,nb);
Yff=zeros(nb,nb);
Yft=zeros(nb,nb);
Yft=zeros(nb,nb);
for k=1:nl   %处理各支�?    
    I=branch(k,1);  
    J=branch(k,2);   
    stat = branch(k, BR_STATUS); 
    Zt=branch(k,BR_R)+j*branch(k,BR_X);          % 
    if Zt~=0
        Yt=stat/Zt;                       % t
    end  
    Ym=j*branch(k,5)/2;  
    Ytt=Yt+Ym;    
    ktap=branch(k,9);       %ratio of transformer
    tap=ktap;      
    if (tap==0)&(J~=0)                   % 普�?支路
        Ybr(I,I)=Ybr(I,I)+Ytt;       
        Ybr(J,J)=Ybr(J,J)+Ytt;       
        Ybr(I,J)=Ybr(I,J)-Yt;       
        Ybr(J,I)=Ybr(I,J);   
    end   
    if (tap==0)&(J==0)                   % 对地支路
        Ybr(I,I)=Ybr(I,I)+Ym;   
    end      
    if tap>0                                  %变比大于0的变压器支路    
        Ybr(I,I)=Ybr(I,I)+Yt/tap/tap;       
        Ybr(J,J)=Ybr(J,J)+Ytt;       
        Ybr(I,J)=Ybr(I,J)-Yt/tap;       
        Ybr(J,I)=Ybr(I,J);   
    end      
    if tap<0                                % 变比小于0的变压器支路      
        Ybr(I,I)=Ybr(I,I)+Ytt;       
        Ybr(J,J)=Ybr(J,J)+tap*tap*Yt;       
        Ybr(I,J)=Ybr(I,J)+tap*Yt;       
        Ybr(J,I)=Y(I,J);   
    end
end
Ybrr=sparse(Ybr);
for busCount=1:nb
    Ygb(busCount,busCount)=(Ygb(busCount,GS)+j*bus(busCount,BS))/baseMVA; %
end
Ygbb=sparse(Ygb);
Ybus=(Ybr+Ygb);%transform Ybus into sparse matrix
%matpower code
f = branch(:, F_BUS);                           %% list of "from" buses
t = branch(:, T_BUS);                           %% list of "to" buses
stat = branch(:, BR_STATUS);                    %% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + j * branch(:, BR_X));   %% series admittance
Bc = stat .* branch(:, BR_B);                           %% line charging susceptance
tap = ones(nl, 1);                              %% default tap ratio = 1
i = find(branch(:, TAP));                       %% indices of non-zero tap ratios
tap(i) = branch(i, TAP);                        %% assign non-zero tap ratios
tap = tap .* exp(-j*pi/180 * branch(:, SHIFT)); %% add phase shifters
Ytt = Ys + j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;
if nargout > 1
    i = [[1:nl]'; [1:nl]'];     %% double set of row indices    
    Yf = sparse(i, [f; t], [Yff; Yft]);
    Yt = sparse(i, [f; t], [Ytf; Ytt]);
end
%%%%

