function [Ra,Rb,Na,Nb,Ta,Tb] = decompose_Homography(H)

I=eye(3);

S=H'*H-I;

M=zeros(3);

M(1,1) = S(2,3)*S(3,2)-S(2,2)*S(3,3);
M(2,1) = S(1,3)*S(3,2)-S(1,2)*S(3,3);
M(3,1) = S(1,3)*S(2,2)-S(1,2)*S(2,3);

M(1,2) = S(2,3)*S(3,1)-S(2,1)*S(3,3);
M(2,2) = S(1,3)*S(3,1)-S(1,1)*S(3,3);
M(3,2) = S(1,3)*S(2,1)-S(1,1)*S(2,3);

M(1,3) = S(2,2)*S(3,1)-S(2,1)*S(3,2);
M(2,3) = S(1,2)*S(3,1)-S(1,1)*S(3,2);
M(3,3) = S(1,2)*S(2,1)-S(1,1)*S(2,2);

check=[M(1,2)^2-M(1,1)*M(2,2);
    M(1,3)^2-M(1,1)*M(3,3);
    M(2,3)^2-M(2,2)*M(3,3)];

fprintf("%f\n",check);

[~,idx]=max([abs(S(1,1)),abs(S(2,2)),abs(S(3,3))]);

Na(3)=NaN;
Nb(3)=NaN;
vs=2*((1+trace(S))^2+1-trace(S^2));
v=sqrt((vs));

Te = norm(sqrt((2+trace(S)-v)));
rho = (sqrt((2+trace(S)+v)));


if idx==1
    Na=[S(1,1);
        S(1,2)+sqrt(M(3,3));
        S(1,3)+sn(M(2,3))*sqrt(M(2,2))];
    Nb=[S(1,1);
        S(1,2)-sqrt(M(3,3));
        S(1,3)-sn(M(2,3))*sqrt(M(2,2))];
    
elseif idx==2
    Na=[
        S(1,2)+sqrt(M(3,3));
        S(2,2);
        S(2,3)-sn(M(1,3))*sqrt(M(1,1))];
    Nb=[
        S(1,2)-sqrt(M(3,3));
        S(2,2);
        S(2,3)+sn(M(1,3))*sqrt(M(1,1))];
elseif idx==3
    Na=[
        S(1,3)+sn(M(1,2))*sqrt(M(2,2));
        S(2,3)+sqrt(M(1,1));
        S(3,3);
        ];
    Nb=[
        S(1,3)-sn(M(1,2))*sqrt(M(2,2));
        S(2,3)-sqrt(M(1,1));
        S(3,3);
        ];
end

Na= Na./norm( Na);
Nb= Nb./norm( Nb);

Tas=(Te/2)*(sn(M(idx,idx))*rho*Nb - Te*Na);
Tbs=(Te/2)*(sn(M(idx,idx))*rho*Na - Te*Nb);

Ra = H*(I-(2/v)*Tas*(Na'));
Rb = H*(I-(2/v)*Tbs*(Nb'));

Ta=Ra*Tas;
Tb=Rb*Tbs;

end

function s=sn(val)
if val>=0
    s=1;
else
    s=-1;
end
end