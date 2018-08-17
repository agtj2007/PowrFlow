   
function [JJ] = JaccobiMatrix_LB110(Va,Vm,nb,npq,npv,G,B,P,Q,pv,pq,ref)
    npvpq=npv+npq;
    H=zeros(npvpq,npvpq);
    N=zeros(npvpq,npq);
    M=zeros(npq,npvpq);
    L=zeros(npq,npvpq);
    count1=1;
    count2=1;
    count3=1;
    count4=1;  
%calculate H matrix and N matrix
    i=[pv' pq']; 
         for ii=i
             for kk=1:nb
                 if ii==kk
                     H(count1,kk)=-1*(Vm(ii)*Vm(kk)*B(ii,kk)+Q(ii)); 
                     N(count1,kk)=Vm(ii)*Vm(kk)*G(ii,kk)+P(ii);  
                 else
                     H(count1,kk)=-1*(Vm(ii)*Vm(kk)*(B(ii,kk)*cos(Va(ii)-Va(kk))-G(ii,kk)*sin(Va(ii)-Va(kk))));        
                     N(count1,kk)=Vm(ii)*Vm(kk)*(G(ii,kk)*cos(Va(ii)-Va(kk))+B(ii,kk)*sin(Va(ii)-Va(kk)));
                 end                                                
              end
          count1=count1+1;
         end 
         H1=H(:,[pv;pq]); 
         NJ=N(:,pq);  
%calculate M matrix and L matrix
     i=pq';
     for ii=i
         for kk=1:nb
             if ii==kk
                 M(count3,kk)=-Vm(ii)*Vm(kk)*G(ii,kk)+P(ii);
                 L(count3,kk)=-Vm(ii)*Vm(kk)*B(ii,kk)+Q(ii);  
             else                
                 M(count3,kk)=-1*(Vm(ii)*Vm(kk)*(G(ii,kk)*cos(Va(ii)-Va(kk))+B(ii,kk)*sin(Va(ii)-Va(kk))));
                 L(count3,kk)=-1*(Vm(ii)*Vm(kk)*(B(ii,kk)*cos(Va(ii)-Va(kk))-G(ii,kk)*sin(Va(ii)-Va(kk))));                 
             end                                  
         end
         count3=count3+1;
     end
     MJ=M(:,[pv;pq]);   
     LJ=L(:,pq);
     JJ=[H1 NJ;
        MJ LJ];   
end

    
    