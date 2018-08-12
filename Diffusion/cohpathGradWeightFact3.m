% Coherence pathway Gradient weighting factor calculation
% calculate the gradient weighting factors bij
% see e.g. Sigmund2006, Song2004
% by Steven Baete, July 2011

% CP    contains coherence pathways, each row one coherence pathway, 0, +1
% or -1
% tau   timing between consecutive pulses, the first time is of course 0,
% the second time is the time between the first and the second pulse
% Grads   array of gradient objects
% intstep integration step, e.g. 1 ms (this should normally be fine enough)
% gamma gyromagnetic ratio of your favorite nucleus 
% B     contains the gradient weighting factors in s/m^2
%       first index: echopulsenumber
%       second index: bxx, byy, bzz, 2*bxy, 2*byz, 2*bxz
%       Sigmund2006 p9 equation (5) and (6) and matrix filling as in (9)
% F     senitivity of the phase shift to the flow in s/m
%       se Cho2007b, JMR 186, p11-16
%       Raguin2007, JMR 184, p337-343
% q contains the magnetisation orientations on the Gt timings (in three
% directions (see output of cohpathMagnOrient2)
%   first index: coherence pathway
%   second index: time (in steps intstep)
% version 2: gradients are passed as objects from the class GradPulse

% CP
% tau
% grads
% intstep
% gamma      [units....]
% q
% silent     false
function [B,F] = cohpathGradWeightFact2(CP,tau,grads,intstep,gamma,q,silent)

try
    a = silent;
catch
    silent = false;
end;

q = single(q);

n=find(CP(1,:)==-100,1,'first')-2;
if (length(tau) ~= n)
    error('Length of timing vector wrong');
end;

%tauinc = cumsum(tau); tauinc(n+1)=tauinc(n)*100;

% construct gradient matrices, these contain gradient values for each
% direction, for each timepoint

maxTE = max(CP(:,n+3)+ 2*tau(2));
nsteps = ceil(maxTE/(intstep/1000));

G = single(zeros(3,nsteps)); % in T/m
for j=1:length(grads)
    if (floor((grads{j}.StartTime + grads{j}.TotalTime)/(intstep*1000)) < nsteps)
        try
            axis = grads{j}.Axis;
            StartT = grads{j}.StartTime;
            TotalT = grads{j}.TotalTime;
            tt = round(StartT/(intstep*1000)):(round((StartT+TotalT)/(intstep*1000))-1);
            G(axis,tt)=G(axis,tt)  + GradShape(grads{j},intstep*1000);
        catch
%             fig2=figure;
%             cohpathPlotGradients(grads,tau*1e6,tau,fig2);
            for jj=1:j
                grads{jj}
            end;
            size(G)
            grads{j}.Axis
            min( round(grads{j}.StartTime/(intstep*1000)):( round(grads{j}.StartTime/(intstep*1000)) + length( grads{j}.GradShape(intstep*1000) ) -1 ) )
            max( round(grads{j}.StartTime/(intstep*1000)):( round(grads{j}.StartTime/(intstep*1000)) + length( grads{j}.GradShape(intstep*1000) ) -1 ) )
            G(grads{j}.Axis,max(1,round(grads{j}.StartTime/(intstep*1000)):( round(grads{j}.StartTime/(intstep*1000)) + length( grads{j}.GradShape(intstep*1000) ) -1 ))) = NaN;
            display('Yup, created a NaN !');
        end;
    end;
end;

%% Trans A >> P
% G = -G;
% G(1,:) = -G(1,:);

%save('/Exp/Meditate/gvectorMatlab','G');

% plot to check G-filling
% fig=figure;hold on;
% % cohpathPlotGradients(grads,tau(2)*1e6,tau,fig);
% for j=1:3
%     plot(intstep*(((1:nsteps)+1/2)/1000),G(j,:)*0.4/0.015/0.535*0.301-0.5-(j-1),'-r');
% end;
% 
% for j=1:3
%     fig=figure;hold on;
%     % cohpathPlotGradients(grads,tau(2)*1e6,tau,fig);
%     plot(intstep*(((1:nsteps)+1/2)/1000),G(j,:),'-r');
% end;
% 
% error

% j=4;
% if( size(q,1) >= j)
%     figure;plot(q(j,:));hold on;
%     plot((((1:nsteps)+1/2)-(tau(2)*1e6)/intstep/1000-1),G(1,:)*1000,'-r');
%     plot((((1:nsteps)+1/2)-(tau(2)*1e6)/intstep/1000-1),G(2,:)*1000,'-g');
%     plot((((1:nsteps)+1/2)-(tau(2)*1e6)/intstep/1000-1),G(3,:)*1000,'-k');
%     ylabel('G [mT/m]');xlabel('t [\mu s]');
%     TE = CP(j,n+3);
%     nsteps = ceil(TE/(intstep/1000));
%     plot(q(j,1:nsteps).*G(3,(1:nsteps)+(tau(2)*1e6)/intstep/1000-1)*1000,'-m');
% end;

% G([2,3],:) = G([3,2],:);

B = zeros(size(CP,1),6);
F = zeros(size(CP,1),3);
    
lq = size(q,1);
lG = size(G,1);
% loop over the echoes
for j=1:size(CP,1)
    TE = CP(j,n+3);
    nsteps = floor(TE/(intstep/1000));
      
    kk = single(zeros(3,nsteps)); 
    tt = 1:nsteps;
    q2 = q(j,tt);
    G2 = G(:,(tt)+(tau(2)*1e6)/intstep/1000-1);
    for k=1:3
% old        temp = q(j,tt).*G(k,(tt)+(tau(2)*1e6)/intstep/1000-1);
%         kk(k,:) = single(cumsum(temp,2)); % in T/m
%         F(j,k) = sum(temp.*(tt));
        temp = q2.*G2(k,:); %faster
        kk(k,:) = single(cumsum(temp,2)); % in T/m
        if (~silent)
            F(j,k) = sum(temp.*(tt));
        end;
    end;

	% [nbe] Save q: instataneous coherence pathway
	%           kk: accumulated integral of the gradient
%     if ( j == size(CP,1) && size(CP,1)==17)
%         save('/Exp/Meditate/qvector2','q');
%         save('/Exp/Meditate/kvectorMatlab','kk');
%     end;
    %kk = gamma * intstep/1000 * kk;  % in T/m * 1/(s*T) *s =  1/m
    
%      figure; hold on;
%     plot(intstep*(1:size(kk,2))/1000,kk(1,:),'-r');
%     plot(intstep*(1:size(kk,2))/1000,kk(2,:),'-b');
%     plot(intstep*(1:size(kk,2))/1000,kk(3,:),'-g');
%     plot(intstep*(1:size(q,2))/1000,q(j,:)*1000,':r');
%     display([' Echo ' num2str(j) ': [' num2str(kk(1,nsteps)) ', ' num2str(kk(2,nsteps)) ', ' num2str(kk(3,nsteps)) ']']);
% 
%      figure; hold on;
%     plot(intstep*(1:size(kk,2))/1000,kk(1,:).*kk(2,:),'-r');
%     plot(intstep*(1:size(kk,2))/1000,kk(2,:).*kk(3,:),'-b');
%     plot(intstep*(1:size(kk,2))/1000,kk(1,:).*kk(3,:),'-g');
%     plot(intstep*(1:size(q,2))/1000,q(j,:)*1000,':r');
    
    % are the directions polarities correct?
%     kk(1,:)=-kk(1,:);
%     kk(2,:)=-kk(2,:);
%     kk(3,:)=-kk(3,:);
    % are the directions correct?
%     kk([1,2],:)=kk([2,1],:);
%     kk([1,3],:)=kk([3,1],:);
%     kk([2,3],:)=kk([3,2],:);

    % xx, yy, zz
    fact = gamma * gamma * intstep/1000 * intstep/1000 * intstep/1000;
    tt = 1:nsteps;
    kk = kk(:,tt);
    for k=1:3
        B(j,k) = fact * sum(kk(k,:).*kk(k,:));        % in (1/m)^2*s
    end;
    B(j,4) = fact * 4*sum(kk(1,:).*kk(2,:));       
    B(j,5) = fact * 4*sum(kk(2,:).*kk(3,:));       
    B(j,6) = fact * 4*sum(kk(1,:).*kk(3,:));   
%     clear kk;
end;

if (~silent)
    fact = gamma*intstep/1000/(nsteps)*TE;
    F = fact*F; % in (rad Hz/T * s *T/m *s: rad s/m )
end;

B = double(B);
