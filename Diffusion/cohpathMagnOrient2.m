% Coherence pathway magnetisation orientation at discrete timesteps
% calculate the diffusion gradient strengths determined by the first three
% gradients (in three directions) and the echo formation condition
% see e.g. Sigmund2006 p9
% by Steven Baete, July 2011

% CP    contains coherence pathways, each row one coherence pathway, 0, +1
% or -1
% tau   time interval between (the center of) consecutive pulses, the first time is of course 0,
% the second time is the time between the first and the second pulse
% q contains the magnetisation orientations on the Gt timings (in three
% directions 
%   first index: coherence pathway
%   second index: time (in steps intstep)

% [nbe]
% tau     : array  : [us] : E.g., for TSE [0, TE, TE, ..., TE]
% intstep : scalar : [us] : E.g., 10us
function q = cohpathMagnOrient2(CP,tau,intstep)

n=find(CP(1,:)==-100,1,'first')-2;
if (length(tau) ~= n)
    error('Length of timing vector wrong');
end;
tauinc = cumsum(tau); tauinc(n+1)=tauinc(n)*100;

maxTE = max(CP(:,n+3));
nsteps = ceil(maxTE/intstep);
q = zeros(size(CP,1),nsteps);
% q2 = zeros(1,nsteps);
    
    
% loop over the echoes
for j=1:size(CP,1)
    TE = CP(j,n+3);
    nsteps = ceil(TE/intstep);
    % fill in the magnetisation orientations at the grid points

%     l=1;  % slow version replaced for a faster one beneath
%     for m=1:length(tauinc)
%         while ((intstep*(l+0.5) < tauinc(m)) && l <= nsteps)
%             q(j,l) = CP(j,m);
%             l = l+1;
%         end;
%     end;
    l=1;
    for m=1:length(tauinc)
        % old and slow
%         lend = min(find( (intstep*((1:nsteps)+0.5)) >= tauinc(m),1,'first' ) +1,nsteps);
        % proposed
        lend = max(min( ceil((tauinc(m)/intstep-0.5)) + 1, nsteps),2);
%         if (lend ~= lend2 )
%             display(['jamojamo, lend = ' num2str(lend) ' and lend2 = ' num2str(lend2)]);
%         end;
        if (isempty(lend))
            lend=nsteps;
        end;
        if ((lend == -1))
            lend=nsteps;
        end;
        q(j,l:lend) =  CP(j,m); %old
%         q2(l:lend) = CP(j,m); %proposed
%        q(j+((l:lend)-1)*size(CP,1)) =  CP(j,m); %proposed and slow
%        q(j,l:lend) =  CP(j,m) * ones(1,lend-l+1); %proposed and slow
%         if ((l > 0) && m > 2 )
%             q(j,(l)) = 0;
%         end;
        l=lend+1;
    end;
%     q(j,:) = q2; %proposed
end;