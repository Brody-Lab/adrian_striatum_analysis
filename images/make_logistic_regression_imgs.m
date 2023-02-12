P=get_parameters;

cov=0.6;
sigma=ones(3)*cov;
sigma = sigma + eye(3)*(1-cov);
vec=[1 1 1];
dist=0.9;
mu1=vec*dist;
mu2=-vec*dist;

n=50;
X=mvnrnd(mu1,sigma,n);
X2=mvnrnd(mu2,sigma,n);

figure('color','w');
h(1)=scatter3(X(:,1),X(:,2),X(:,3),20,'k','filled');hold on;
h(2)=scatter3(X2(:,1),X2(:,2),X2(:,3),20,'k');
line_length=4;

quiver3(0,0,0,vec(1)*line_length,vec(2)*line_length,vec(3)*line_length,'LineWidth',2,'Color','k','MaxHeadSize',0.6);
text(4.5,4.5,4.5,'$\vec{\beta}$','FontSize',15,'Interpreter','latex');
text(X(1,1),X(1,2),X(1,3),'$R(t)$','Interpreter','latex','color','r');

b_hat = vec./norm(vec);
proj = sum(X(1,:).*b_hat).*b_hat;
rej = X(1,:)-proj;

line([X(1,1) rej(1)],[X(1,2) rej(2)],[X(1,3) rej(3)],'color','r','LineStyle',':','LineWidth',2);

text(mean([X(1,1) rej(1)]),mean([X(1,2) rej(2)]),mean([X(1,3) rej(3)]),'$DV(t)$','Interpreter','latex','color','r');


sz=4;
patch_pos=null(vec);
patch_pos = [ [0 0 0]' patch_pos sum(patch_pos,2)];
patch_pos = patch_pos - mean(patch_pos,2);
patch_pos = patch_pos *sz;
patch_pos = patch_pos(:,[1 2 4 3]);



patch('XData',patch_pos(1,:),'YData',patch_pos(2,:),'ZData',patch_pos(3,:),'FaceAlpha',0.5);
lim=[-4 4];
set(gca,'xlim',lim,'ylim',lim,'zlim',lim);
axis square
xlabel('R_1');
ylabel('R_2');
zlabel('R_3');
set(gca,'xtick',[],'ytick',[],'ztick',[],'box','off',P.axes_properties{:});

scatter3(X(1,1),X(1,2),X(1,3),30,'r','filled');hold on;
g=legend(h,{'R choice','L choice'},'Location','Northeast');



mdl =  fitglm([X;X2],[ones(n,1);zeros(n,1)],'linear','Distribution','binomial');

DV=mdl.Fitted.LinearPredictor;


figure('color','w');
xs=linspace(-10,10,100);
plot(xs,(1./(1+exp(-xs))),'LineWidth',2,'color','k');hold on;
scatter(DV(1:n),ones(n,1),20,'k','filled');
scatter(DV(1),1,30,'r','filled');
scatter(DV(n+1:end),zeros(n,1),20,'k');
set(gca,'xlim',[-5 5],'box','off',P.axes_properties{:});
xlabel('$DV(t)$','Interpreter','latex');
ylabel('$p(R)$','Interpreter','latex');


