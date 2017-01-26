
function ee = comparch2(type,L,n,iter)

if nargin < 4, iter = 25; end
if nargin < 3 | isempty(n), n = 50; end
if nargin < 2, L = []; end

W = classfs;
if ~isempty(L), W = W(L); end

rand('state',0);
ee.names = getname(W{1});
for j=2:length(W)
	ee.names = char(ee.names,getname(W{j}));
end

learnsizes = n;
	
e = zeros(length(W),length(learnsizes),iter);
b = genarche(type,1000,0);
for j=1:iter
	for m = 1:length(learnsizes)
		a = genarche(type,learnsizes(m),j);
		for i=1:length(W)
			e(i,m,j) = b*(a*W{i})*testc;
			disp([j,m,i])
		end
	end
end

ee.error = e;
ee.iter = iter;
ee.learnsizes = learnsizes;

disp(type)
disp(' ')
for m = 1:length(learnsizes)
	em = mean(squeeze(e(:,m,:))');
	es = std(squeeze(e(:,m,:))');
	[eem,J] = sort(em);
	for i=1:length(W)
		fprintf(1,'%7.4f %8.5f   %s \n',em(J(i)),es(J(i))/sqrt(iter),ee.names(J(i),:));
	end
end

if isempty(L) % full set, store result
	save([type '_' num2str(iter) '_' num2str(learnsizes(1))],'ee');
end

%figure;
%scatterd(genarche(type,200,0));
%axis equal
%fsave('scatter_plot',14,2,5,1.5);
%ee.error = mean(e,3);
%ee.xvalues = learnsizes;
%ee.std = std(e,[],3)/sqrt(iter);
%plotr(ee,[],[],[],'errorbar')
%save(['arch2_' getmapping_file(W{n})],'ee','iter');

			
