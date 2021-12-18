%%
% Create the connection matrix used in the diffusion approximation such
% that the in-degrees are fixed
% Note: There are no feedforward connections, since these are approximated
% in the main function
%%
% function [wind_pre,wipre,wstr_pre,pc_pre] = Diff_approx_gen_weights_v2(Ncells,p0,J,pinds)
function [wind_post,wipost,wstr_post] = Diff_approx_gen_weights_v2(Ncells,p0,J,pinds)

%% set up recurrent weight matrix
Ntot = sum(Ncells);
Npop = length(Ncells);

Maxw = round(Ntot*Ntot*0.3); % maximum number of weights in the weight matrix
wind_pre = zeros(Ntot+1,1); % column of w corresponding to the start of the ith neuron's projections                
	
wipre = zeros(Maxw,1);
wstr_pre = zeros(Maxw,1);

% pc_pre = zeros(Maxw,1);
pc_vec = zeros(Ntot,1);
for kk = 1:Npop
    pc_vec(pinds(kk):pinds(kk+1)-1) = kk;
end

syncount = 1;

for cc = 1:Ntot
    
    pc = pc_vec(cc);
    
    wind_pre(cc) = syncount;
    
    for pp = 1:Npop
    
        % probability of a connection
        prob = p0(pp,pc);

        % Method 2: fix the total number of connections to be Ncells(pp)*prob
        iconns = randperm(Ncells(pp),round(Ncells(pp)*prob))+sum(Ncells(1:(pp-1)));
        
        wipre(syncount:(syncount+length(iconns)-1)) = iconns;
        
        
        wstr_pre(syncount:(syncount+length(iconns)-1)) = J(pp,pc);
        
%         pc_pre(syncount:(syncount+length(iconns)-1)) = pp;
            
        syncount = syncount+length(iconns);
    end
     
end

wind_pre(Ntot+1)= syncount-1;
% reshape these vectors to get rid of unnecessary zeros
wipre = wipre(1:(syncount-1));
wstr_pre = wstr_pre(1:(syncount-1));
% pc_pre = pc_pre(1:(syncount-1));



%% Now turn it around
wipost = zeros(size(wipre));
for cc = 1:Ntot
    wipost(wind_pre(cc):wind_pre(cc+1)) = cc;
end
[~, new_order]= sort(wipre);

wipost = wipost(new_order);
wind_post = find(diff(wipost)<0)+1;
wind_post = [1; wind_post; size(wipost,1)];

wstr_post = wstr_pre(new_order);


% syncount = 1;
% % loop through the populations
% for pp = 1:Npop
%     
%     % loop through each neuron
%     for cc=1:Ncells(pp)
%         
%         starting_Index = sum(Ncells(1:(pp-1)));
%         wind(cc + starting_Index) = syncount;
%         
%         % find which neurons are connected to neuron cc
%         for qq = 1:Npop
%         
%             % probability of a connection
%             prob = p0(pp,qq);
%             
%             % Method 1: flip weighted coins!            
% %             iconns = find(rand(Ncells(qq),1) < prob) + sum(Ncells(1:(qq-1)));
%             
%             % Method 2: fix the total number of connections to be Ncells(qq)*prob
%             iconns = randperm(Ncells(qq),round(Ncells(qq)*prob))+sum(Ncells(1:(qq-1)));
%              
%             % record the connections
%             wipre(syncount:(syncount+length(iconns)-1)) = iconns;
%             wstr(syncount:(syncount+length(iconns)-1)) = J(pp,qq);
%             
%             syncount = syncount+length(iconns);
%         end
%     end
% end



end