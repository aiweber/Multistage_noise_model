
function likes = mult_likelihood_wrapper(g,upNoise,multNoise,ys,nlPrs,lL,lU,yU,lCut,flag)
% likes = mult_likelihood_wrapper(g,upNoise,multNoise,ys,nlPrs,lL,lU,yU,lCut,flag)
%
% Wrapper for mult_likelihood that enables it to operate on vector inputs

likes = zeros(size(ys));
for yValIdx = 1:length(ys)
    y = ys(yValIdx);
   
    %%%% adjust integration bounds for each y
    if flag
        lL2 = max(lCut,y-5*sqrt(yU)*multNoise);
        lU2 = min(yU,lL2+10*sqrt(yU)*multNoise);
    else
        lL2 = lL;
        lU2 = lU;
    end
    
    likes(yValIdx) = mult_likelihood([g upNoise multNoise y nlPrs],lL2,lU2);
    
end
