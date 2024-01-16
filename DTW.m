clear
close all
clc

tic
%% Import data
% seismic data
load('pp.mat');   % P-waves
load('ps.mat');   % S-waves

%% S-waves preprocessing
for m=1:size(ps,2)
    for n=1:(size(ps,1)-1)/2
        if abs(ps(2*n-1,m))>abs(ps(2*n,m))
            ps1(n,m)=ps(2*n-1,m);
        else
            ps1(n,m)=ps(2*n,m);
        end
    end
end
ps1(2501,:)=ps(5001,:); % ps1 is the preprocessed S-waves

%% Registration
x=length(pp);
y=length(ps1);
L=36;   % Slope constraint
for kkk=1:1     % Number of iterations
    for k=1:1601   % Trace number
        %% Alignment errors
        for i=1:x   % Set to infinity if not within path constraints
            for j=1:y
                if abs(i-j)>L
                    e(i,j)=Inf;
                    continue;
                end
                e(i,j)=(pp(i,k)-ps1(j,k))^2;
            end
        end
        
       %% Accumulation
        d(1,:)=e(1,:);
        for i=2:x       % i represents the moment of the P-waves
            for j=1:y   % j represents the moment of the S-waves
                if abs(i-j)>L
                    d(i,j)=Inf;
                    continue;
                end
                % d(i,j)=e(i,j)+min(d(i-1,j-1),d(i-1,j),d(i-1,j+1));
                if j==1
                    d(i,j)=e(i,j)+d(i-1,j);
                    continue;
                end
                if j==2
                    d(i,j)=e(i,j)+min(d(i-1,j-1),d(i-1,j));
                    continue;
                end
                a=min([d(i-1,j-2) d(i-1,j-1) d(i-1,j)]);
                d(i,j)=e(i,j)+a;
            end
        end
        
       %% Backtracking
        u(i,k)=0;
        for j=1:x-L-2
            i=i-1;
            a=min([d(i,i+u(i+1,k)-1) d(i,i+u(i+1,k)) d(i,i+u(i+1,k)+1)]);
            if a==d(i,i+u(i+1,k))           
                u(i,k)=u(i+1,k);
            elseif a==d(i,i+u(i+1,k)+1)
                u(i,k)=u(i+1,k)+1;
            elseif a==d(i,i+u(i+1,k)-1)
                u(i,k)=u(i+1,k)-1;
            end
        end

        for i=1:x
            s(i,k)=ps1(i+u(i,k),k);
        end
        c2(k)=sum(pp(:,k).*s(:,k))/sqrt(sum(pp(:,k).^2)*sum(s(:,k).^2));   % Similarity
    end
end
timecost=toc


%% Figure
% Comparison diagram of single trace
% k=1;    % Number of trace
% subplot(311);
% plot(pp(:,k),'LineWidth', 1.2);
% t1=title('(a)');
% t3.FontSize=24;
% set(legend,...
%     'Position',[0.766840273304325 0.85309882770581 0.137890627495945 0.0457286441865278]);
% xlim([0,2501])
% ylim([-60000,80000])
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Amplitude')
% xlabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')
% 
% subplot(312);
% plot(ps1(:,k),'LineWidth', 1.2);
% t2=title('(b)');
% t3.FontSize=24;
% set(legend,...
%     'Position',[0.767361106637658 0.552596315142996 0.137890627495944 0.0457286441865279]);
% xlim([0,2501])
% ylim([-60000,80000])
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Amplitude')
% xlabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')

% subplot(313);
% plot(s(:,k),'LineWidth', 1.2);
% t3=title('(c)');
% t3.FontSize=24;
% set(legend,...
%     'Position',[0.675434022205363 0.2541038528554 0.230078129693866 0.0437185938873482]);
% xlim([0,2501])
% ylim([-60000,80000])
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Amplitude')
% xlabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')

%% Full trace diagram
% figure 
% imagesc(s)
% colormap(gray)
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')
% xlabel('\fontname{Times New Roman}\fontsize{24}Trace')

% figure
% imagesc(pp)
% colormap(gray)
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')
% xlabel('\fontname{Times New Roman}\fontsize{24}Trace')

% figure
% imagesc(ps)
% colormap(gray)
% set(gca,'FontSize',24,'fontname','Times New Roman')
% ylabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')
% xlabel('\fontname{Times New Roman}\fontsize{24}Trace')

%% Offset map  u=j-i
% k=1
% u1=u(:,k);
% pp1=pp(:,k);
% ps11=ps1(:,k);    
% pp11=pp1+140000;    
% figure('Name','Offset map')
% plot(pp11,'k')
% hold on
% plot(ps11,'k')
% xlim([0,2501])
% xlim([0,x])
% ylim([-60000,180000])
% q=5;
% for i=1:floor(x/5)
%     a=[q*i,q*i+u1(q*i)];
%     b=[pp11(q*i),ps11(q*i+u1(q*i))];
%     plot(a,b,'k')
% end
% set(gca,'FontSize',24,'fontname','Times New Roman')
% xlabel('\fontname{Times New Roman}\fontsize{24}Time(ms)')

% hold off





