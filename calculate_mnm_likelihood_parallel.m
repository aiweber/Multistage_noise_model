function [loglike,loglikes] = calculate_mnm_likelihood_parallel(prs,gSub,rSub,fname)
% [loglike,loglikes] = calculate_mnm_likelihood_parallel(prs,gSub,rSub,fname)

disp(num2str(prs))
upNoise = prs(1);   % std of upstream noise
multNoise = prs(2); % strength of multiplicative noise (scales with output of nonlinearity)
downNoise = prs(3); % std of downstream noise
nlPrs = prs(4:7);
downP = prs(8); % probability that downstream noise is drawn from Gaussian distribution (noise is 0 with probability 1-downP)

fid = fopen([fname '.txt'],'r');
a = fscanf(fid,'%f %f %f %f %f %f %f %f',[8 Inf]);
fclose(fid);

a = [a prs'];

fid = fopen([fname '.txt'],'w');
fprintf(fid,'%f %f %f %f %f %f %f %f\n',a);
fclose(fid);

nStds = 4; % number of standard deviations of each distribution to capture
derivCut = 1e-3; % cutoff for not calculating integral (where nonlinearity is very flat)
d = nlPrs(4);

if nlPrs(1)*nlPrs(2)<= derivCut % if nonlinearity if flat
    gCut = max(gSub)+nStds*upNoise; % all prob lies in flat region
else
    gCut = (log(nlPrs(1)*nlPrs(2)/derivCut -1 )+nlPrs(3))/-nlPrs(2);
end
lCut = nl_sr(nlPrs,gCut);

%%% checks that nonlinearity/noise params are reasonable 
if nl_sr(nlPrs,max(gSub)+5*upNoise) == Inf % if largest output of nonlinearity is Inf
    loglike = -Inf;
    loglikes = -Inf*ones(size(gSub));
    return
end

% make sure largest response is possible at nl_sr(g), given upstream noise
g = max(gSub);
rObs = max(rSub); % not corresponding points, but make sure max r possible near max g
maxPossR = nl_sr(nlPrs,g+nStds*upNoise) + nStds*multNoise*sqrt(nl_sr(nlPrs,g+nStds*upNoise)) + nStds*downNoise;
maxPossR = max(maxPossR,0);  % can't be negative
if rObs>maxPossR
    loglike = -Inf;
    loglikes = -Inf*ones(size(gSub));
    return
end


%%
%%%%%%%%%
% assess likelihood for each point rObs given g
%%%%%%%%%

loglikes = zeros(length(gSub),1);
parfor i = 1:length(gSub)
    g = gSub(i);
    rObs = rSub(i);

    if g+nStds*upNoise <= gCut  %%% if all probability in flat region of NL, output of NL is delta function at nlPrs(4)
        if rObs == 0
            prObs = downP*.5*(1+erf((.5-d)/(sqrt(multNoise^2*d+downNoise^2)*sqrt(2))));  % with down noise, *downP
            prObs = prObs + (1-downP)*.5*(1+erf((.5-d)/(sqrt(multNoise^2*d)*sqrt(2))));  % add in distr w/o down noise, *(1-downP)

        else % rObs is integer, so integrate between rObs-.5 and rObs+.5
            prObs = downP*.5* (erf((rObs+.5-d)/(sqrt(multNoise^2*d+downNoise^2)*sqrt(2))) - erf((rObs-.5-d)/(sqrt(multNoise^2*d+downNoise^2)*sqrt(2))));
            prObs = prObs + (1-downP)*.5* (erf((rObs+.5-d)/(sqrt(multNoise^2*d)*sqrt(2))) - erf((rObs-.5-d)/(sqrt(multNoise^2*d)*sqrt(2))));
        end
        
    else %%% split into flat and non-flat regions
        
        %%% flat region
        % first calculate fraction of total inputs getting cut off by gCut
        fracCut = .5*(1+erf((gCut-g)/(upNoise*sqrt(2))));
        
        %%% non-flat region
        lL = nl_sr(nlPrs,max(gCut,g-nStds*upNoise));
        lU = nl_sr(nlPrs,g+nStds*upNoise);
                
        yL = nStds^2*multNoise^2/4 - 4*nStds*multNoise^2/2;
        yU = lU + nStds*multNoise*sqrt(lU); % this vec will be used once downstram noise is added in too

        if (rObs~=0 && rObs-.5>yU+nStds*downNoise) || rObs+.5<yL-nStds*downNoise
            prObs = 0;
        else
            
            if lU-lL > 10*multNoise
                flag = 1;
            else
                flag = 0;
            end
            
            f = @(ys) mult_likelihood_wrapper(g,upNoise,multNoise,ys,nlPrs,lL,lU,yU,lCut,flag);

            chebfunpref.setDefaults('chebfuneps',1e-4);
            cf = chebfun(f,[yL,yU],'splitting','on');

            nptsy = 10000;
            diffs = (yU-yL)/(nptsy-1);
            
            if fracCut>1e-6
                stdFracCut = multNoise*sqrt(nlPrs(4));
                if stdFracCut/30 > 0 % make sure this isn't zero
                    diffs = min(diffs,max(stdFracCut/30,diffs/100)); % make sure diff is at least as small as 1/30th of std of noise from fraction that was cut (passed through NL at nlPrs(4)),
                                                                     % but make sure it's not so small that yVec becomes huge
                end
                yVec = yL:diffs:yU;
                pyVec = cf(yVec);
                
                % add back in flat region
                pyVec = pyVec/trapz(yVec,pyVec)*(1-fracCut);
                
                if stdFracCut == 0
                    idx = find(abs(yVec-nlPrs(4))==min(abs(yVec-nlPrs(4))));
                    idx = idx(1);
                    cutDist = zeros(size(yVec));
                    cutDist(idx) = fracCut;
                else
                    cutDist = gauss(yVec,nlPrs(4),stdFracCut)/(stdFracCut*sqrt(2*pi));
                    cutDist = cutDist/trapz(yVec,cutDist)*fracCut;
                end
                
                pyVec = pyVec+cutDist;
            else
                yVec = yL:diffs:yU;
                pyVec = cf(yVec);
            end
            pyVec = pyVec/trapz(yVec,pyVec); % normalize to integrate to 1
            cf = [];
            
            %%% pad with 0s for +/- 4 SDs of downstream noise
            yVecAddLow = fliplr(yVec(1):-diffs:yVec(1)-nStds*downNoise);
            yVecAddHigh = yVec(end):diffs:yVec(end)+nStds*downNoise;
            yVec = [yVecAddLow(1:end-1) yVec yVecAddHigh(2:end)];
            pyVec = [zeros(1,length(yVecAddLow)-1) pyVec zeros(1,length(yVecAddHigh)-1)];
            
            %%% convolve with downstream noise to get r (output with all noise, without rounding)
            x = linspace(-(diffs*(length(yVec)-1))/2,(diffs*(length(yVec)-1))/2,length(yVec));
            dn = gauss(x,0,downNoise);
            dn = dn/sum(dn)*downP; % include probability of downstream noise coming from gaussian distribution
            
            Y = fft([pyVec zeros(1,length(dn)-1)]);
            DN = fft([dn zeros(1,length(pyVec)-1)]);
            R = Y.*DN;
            r = ifft(R);

            keepIdx = round(length(r)/2-length(pyVec)/2);
            keepIdx = keepIdx+1:keepIdx+length(pyVec);
            r = r(keepIdx) + pyVec*(1-downP); % add back in component with no noise with probability (1-p)
                                              
            if rObs ~= 0 && rObs-.5>yVec(1) % make sure lower end of rObs-.5 is still within range
                firstV = rObs-.5;
                firstR = interp1(yVec,r,firstV); 
            else
                firstV = yVec(1);
                firstR = r(1);
            end
            if rObs+.5<yVec(end)
                lastV = rObs+.5;
                lastR = interp1(yVec,r,lastV);
            else
                lastV = yVec(end);
                lastR = r(end);
            end
            
            idx = yVec>firstV & yVec<lastV;
            
            vec = [firstV yVec(idx) lastV];
            rKeep = [firstR r(idx) lastR];
            
            prObs = trapz(vec,rKeep);
            
        end
        
    end
    
    loglikes(i) = log(prObs);
end % over all observed g


loglike = sum(loglikes);


