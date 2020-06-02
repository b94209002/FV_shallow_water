function [LL RL LH RH] = BDS_bilinear_constraint(sh,LL,RL,LH,RH)

srt = circshift(sh,[-1,-1]);
srb = circshift(sh,[-1,1]);
sr  = circshift(sh,[-1,0]);
sl  = circshift(sh,[1,0]);
slt = circshift(sh,[1,-1]);
slb = circshift(sh,[1,1]);
st  = circshift(sh,[0,-1]);
sb  = circshift(sh,[0,1]);

minLL = min(min(sb,slb),min(sl,sh));
maxLL = max(max(sb,slb),max(sl,sh));
minRL = min(min(sb,srb),min(sr,sh));
maxRL = max(max(sb,srb),max(sr,sh));
minLH = min(min(st,slt),min(sl,sh));
maxLH = max(max(st,slt),max(sl,sh));
minRH = min(min(st,srt),min(sr,sh));
maxRH = max(max(st,srt),max(sr,sh));

epsilon = 1.e-12;

LL = max(min(LL,maxLL),minLL);
RL = max(min(RL,maxRL),minRL);
LH = max(min(LH,maxLH),minLH);
RH = max(min(RH,maxRH),minRH);

m = size(sh);

for j = 1:m(2)
    for i = 1:m(1)
        temp= [LL(i,j) LH(i,j) RL(i,j) RH(i,j)];
        for k = 1:3
            sumdif = sum(temp) - 4*sh(i,j);
            sgn = sign(sumdif);
            diff  = sgn * (temp - sh(i,j));
            kdp = sum(diff>epsilon);
            if (diff(1) > epsilon)
                [temp(1) kdp sumdif] = BDS_bilinear_iter(temp(1),kdp,sgn,sumdif,minLL(i,j),maxLL(i,j));
            end
            if (diff(2) > epsilon)
                [temp(2) kdp sumdif] = BDS_bilinear_iter(temp(2),kdp,sgn,sumdif,minLH(i,j),maxLH(i,j));
            end
            if (diff(3) > epsilon)
                [temp(3) kdp sumdif] = BDS_bilinear_iter(temp(3),kdp,sgn,sumdif,minRL(i,j),maxRL(i,j));
            end
            if (diff(4) > epsilon)
                [temp(4) kdp sumdif] = BDS_bilinear_iter(temp(4),kdp,sgn,sumdif,minRH(i,j),maxRH(i,j));
            end
        end
        LL(i,j)=temp(1);
        LH(i,j)=temp(2);
        RL(i,j)=temp(3);
        RH(i,j)=temp(4);
    end
end


