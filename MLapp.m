function OUT=MLapp(Z,c,beta1,beta2)
%An all in one call to the 2-parameter MLF, that incorporates the Poincare
%and Power series approximations ONLY to be used when Z>=0 and the MLF is
%one parameter or is two parameter with beta1=beta2
%Z is the input taken and is non-negative
%c is time scale
%beta1 and beta2 are the first and second parameters of the Mittag-Leffler
%function
%This is only to be used for arguments of the form -(cT)^beta for T>0.
%Otherwise the approximations will break down.
if (beta1~=beta2)&&beta2~=1
    error("Need beta1=beta2 or beta2=1")
end
if any(Z<0)
    error("All entries of Z must be non-negative")

end
OUT=zeros(size(Z));
cZb=(c.*Z).^beta1;
Idx1=cZb>=10; %Events to use Poin Asym App.
OUT(Idx1)=AppML(13,beta1,-cZb(Idx1),beta2);

if beta1>=0.5

    Idx2=cZb<=4;

    OUT(Idx2)=PSML(120,beta1,-cZb(Idx2),beta2);

    Idx3=~(Idx1|Idx2);

    OUT(Idx3)=ml(-cZb(Idx3),beta1,beta2,1);

elseif beta1<0.5 && beta1>=0.25

    Idx2=cZb<=2;

    OUT(Idx2)=PSML(200,beta1,-cZb(Idx2),beta2);

    Idx3=~(Idx1|Idx2);

    OUT(Idx3)=ml(-cZb(Idx3),beta1,beta2,1);

    
else

    Idx2=cZb<=1;

    OUT(Idx2)=PSML(100,beta1,-cZb(Idx2),beta2);

    Idx3=~(Idx1|Idx2);

    OUT(Idx3)=ml(-cZb(Idx3),beta1,beta2,1);

   
end

end

function OUT=AppML(M,beta1,Z,beta2)

% OUT=0;
% for II=1:M
% OUT=OUT-(Z.^(-II))./gamma(beta-II.*beta);
% end
%I think below is vectorised version of above

if ~isempty(Z)
OUT=-(gamma(beta2-beta1.*(1:M)).^(-1))*(Z.^(-(1:M).'));
else
    OUT=[];
end

end

function OUT2=PSML(M,beta1,Z,beta2)
% OUT=0;
% for k=0:M
% OUT=OUT+(Z.^k)./gamma(beta2+k.*beta1);
% end
%I think below is vectorised version of above
if ~isempty(Z)
OUT2=real(sum(exp((0:M).*log(Z).'-gammaln(beta2+beta1.*(0:M))),2)).';
OUT2(abs(Z)<=1e-15)=1./gamma(beta2);
else
    OUT2=[];
end
end