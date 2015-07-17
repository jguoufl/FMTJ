function [Ipos Isz]=func_current(energy, HD, AUD, ALD, mu1, mu2,const)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigx=[0 1; 1 0]; sigy=[0 -1i; 1i 0]; sigz=[1 0; 0 -1];
eta=1e-15*i;
ep = energy+eta;
Np=length(HD);

f1=1/(1+exp((energy-mu1)/const.kBT));
f2=1/(1+exp((energy-mu2)/const.kBT));

AD = cell(1,Np);

%% set the blocks for AG=I
for ii = 1:Np
    AD{ii} = ep*eye(2)-HD{ii};
end

%%%% semi-infinite source contact
sig_s=sig_sr2(HD{1},-AUD{1}',eye(2),zeros(2,2),energy);
sig_d=sig_sr2(HD{Np},-AUD{Np-1},eye(2),zeros(2,2),energy);

gam1=1i*(sig_s-sig_s');
gam2=1i*(sig_d-sig_d');


AD{1}  = AD{1} - sig_s;        % add the source self-energy
AD{Np} = AD{Np}- sig_d;       % add the drain  self-energy
% %% computing transmission and current only
% [Grl Grd Gru]=recursealg_concise(Np,ALD,AD,AUD);
% Tr=real(trace((1i*(Grd{Np}-Grd{Np}')-Grd{Np}*gam2*Grd{Np}')*gam2)); %transmission per spin, method 1
% Ipos=Tr*(f1-f2);  % current spectrum
% %%% end of computing transmission and current only

%%% Computing the position-resolved current
con_in=cell(1,Np);
con_out=cell(1,Np);
for ii=2:Np-1  % ballistic
    con_in{ii}=sparse(2,2);
    con_out{ii}=sparse(2,2);
end
con_in{1}=sparse(gam1)*f1;       
con_in{Np}=sparse(gam2)*f2;  
con_out{1}=sparse(gam1)*(1-f1);  
con_out{Np}=sparse(gam2)*(1-f2);
[Grd,Gnl,Gnd, Gpd] = recursealg(Np,ALD,AD,AUD,con_in,con_out);

%[Grl,Grd,Gru,Gnl,Gnd,Gnu,Gpl,Gpd,Gpu,grL,ginL] = recursealg_complete_version(Np,ALD,AD,AUD,con_in,con_out);

for ii=1:Np-1   % position resolved current
    Ipos(ii,1)=trace(imag(AUD{ii}*Gnl{ii}-Gnl{ii}'*ALD{ii}));
    %Ipos(ii)=2*trace(imag(AUD{ii}*Gnl{ii})); % checked same as above
    Ipos(ii,2)=trace(imag(sigz*(AUD{ii}*Gnl{ii}-Gnl{ii}'*ALD{ii})));
    Ipos(ii,3)=trace(imag(sigx*(AUD{ii}*Gnl{ii}-Gnl{ii}'*ALD{ii})));
    Ipos(ii,4)=trace(imag(sigy*(AUD{ii}*Gnl{ii}-Gnl{ii}'*ALD{ii}))); 
end
%%% end of computing the position-resolved current
ii_sz=1; % the position at which current Isz retrieved 
Isz=diag(imag(AUD{ii_sz}*Gnl{ii_sz}-Gnl{ii_sz}'*ALD{ii_sz}));






