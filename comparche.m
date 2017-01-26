
W = {qdc,udc,ldc,nmc,treec,lmnc,parzenc,knnc,knnc([],1),loglc,svc,naivebc,rbnc,fdsc};

ee.names = getname(W{1});
for j=2:length(W)
	ee.names = char(ee.names,getname(W{j}));
end

learnsizes = [2,5,7,10,15,20,30,50,70,100];

iter = 25;
for n = 1:8
	
e = zeros(length(W),length(learnsizes),iter);
b = genarche(W{n},1000,0);
for j=1:iter
	for m = 1:length(learnsizes)
		a = genarche(W{n},learnsizes(m),j);
		for i=1:length(W)
			e(i,m,j) = b*(a*W{i})*testc;
			disp([n,j,m,i])
		end
	end
end

ee.error = mean(e,3);
ee.xvalues = learnsizes;
ee.std = std(e,[],3)/sqrt(iter);
plotr(ee,[],[],[],'errorbar')
save(['arche_' getmapping_file(W{n})],'ee','iter');

end
			
