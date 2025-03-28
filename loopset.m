function L=loopset(input,names)

L.nms=names;

L.adj=sign(input);

E=input;
loopcount = 0;

% The numbers of stages/age-groups, of transitions between them, and the size of independent loop set are reported 
Vertex=size(E,2);               % Number of nodes from the dimension of the transition matrix
Edges=nnz(E);                   % Number of edges from the nonzero elements of the transition matrix
Nullity = Edges-Vertex+1;       % # of independent loops that can be identified. This is also the 
                                %upper limit for # of ecologically-meaningful
                                %loops that can be identified in a given life-cycle graph
fprintf('    Nodes(n)=    %d\n', Vertex)    % display the number of stages/age-groups
fprintf('    Edges(b)=    %d\n', Edges)     % display the number of transitions
fprintf('    Nullity(l)=  %d\n', Nullity)   % display the max number of ecologically-meaningful 
                                            %loops that can be identified

% Identify self-loops, store their elements and elasticity values in fbl and els sets, respectively
fprintf('Identifying the self-loops...\n');
for cntr1=1:size(E,1)
    if (E(cntr1,cntr1) > 0)
        loopcount = loopcount + 1;
        fbl(loopcount,1) = cntr1;
        els(loopcount) = E(cntr1,cntr1);
        E(cntr1,cntr1) = 0;
        fprintf('self-loop %d...\n', loopcount);
    end
end

% Identify two-transition loops, store their elements and elasticity values in fbl and els sets, respectively
fprintf('Identifying the 2-transition loops...\n');
for cntr1=1:size(E,1)
    for cntr2=1:size(E,2)
        if (E(cntr1,cntr2) > 0)&&(E(cntr2,cntr1) > 0)
            loopcount = loopcount + 1;
            fbl(loopcount,1) = cntr1;
            fbl(loopcount,2) = cntr2;
            fprintf('2-transition loop %d - %d\n', cntr1, cntr2);
            if (E(cntr1,cntr2) >= E(cntr2,cntr1))
                E(cntr1,cntr2) = E(cntr1,cntr2) - E(cntr2,cntr1);
                els(loopcount) = E(cntr2,cntr1);
                E(cntr2,cntr1) = 0;
            else
                E(cntr2,cntr1) = E(cntr2,cntr1) - E(cntr1,cntr2);
                els(loopcount) = E(cntr1,cntr2);
                E(cntr1,cntr2) = 0;
            end
        end
    end
end

for bg=3:size(E,1)       % for loop sizes 3 to the size of the matrix
    fprintf('Identifying the %d-transition loops...\n', bg);
%     if bg == size(E,1)
%         for sik=1:bg
%             for kis=1:bg
%                 fprintf('%f\t',E(sik,kis));
%             end
%             fprintf('\n');
%         end
%     end
    et=E';                  % Take transpose of the input matrix-minus-previously-identified 
                            % loops -the function 'reach' below regards rows and columns 
                            % as 'from' and 'to', respectively. 
                            % This is NOT required for system dynamics (dynamic feedback) models!
    a = sign(et);           % Get the adjacency matrix from the projection/elasticity matrix

    R=reach(a);                                 % Calculate the Reachability Matrix
    cyc=cycle(R);

    for i=1:size(cyc,1)							% for each cycle identified
        m=1; y=[];								% initialize loop counter and looplist
        index=find(cyc(i,:))';					% create index-row w/cycle elements
        c_elm=find(cyc(i,:));
        fc=loops(a(c_elm,c_elm),bg);			% id the fb loops of size sz in cycle
	
        % The following for-loop translates the fc back to the original node number
        for j=1:size(fc,1)						% for each fb loop identified
            for k=1:size(fc,2)					% for each element id in the loop
                if (fc(j,k)~=0)					% if there is an element active
                    y(m,k)=index(fc(j,k));	    % extract its coordinates form index
                else							% else, end of fbl ... move to next
                    break
                end
            end
            m=m+1;
        end

%         if ~isempty(y)
%             Edges=nnz(a(c_elm,c_elm));         %Number of edges from the nonzero elements of the transition matrix
%             Vertex=size(c_elm,2);              %Number of nodes from the dimension of the transition matrix
%             Nullity = Edges-Vertex+1;
%             y=ils(y,a,1); 
%             fprintf('    Nodes(n)=    %d\n', Vertex)
%             fprintf('    Edges(b)=    %d\n', Edges)
%             fprintf('    Nullity(l)=  %d\n', Nullity)
%         end;
        
        for cntr1=1:size(y,1)
            min_elast = 999;
            for ig=1:size(y,2)-2       % The lines 166-175 identify the transition 
                if (ig == 1)            % with the smallest elasticity to min_elast
                    dummy_elast = min(E(y(cntr1,ig+1),y(cntr1,ig)), E(y(cntr1,ig+2),y(cntr1,ig+1)));
                else
                    dummy_elast = min(dummy_elast, E(y(cntr1,ig+2),y(cntr1,ig+1)));
                end
                if (dummy_elast < min_elast)
                    min_elast = dummy_elast;
                end
            end
            min_elast = min(E(y(cntr1,1),y(cntr1,size(y,2))), min_elast);
            %         if (E(y(cntr1,1),y(cntr1,2)) > E(y(cntr1,2),y(cntr1,3)))
            %             min_elast = E(y(cntr1,2),y(cntr1,3));
            %         else
            %             min_elast = E(y(cntr1,1),y(cntr1,2));
            %         end
            %         if (min_elast > E(y(cntr1,3),y(cntr1,1)))
            %             min_elast = E(y(cntr1,3),y(cntr1,1));
            %         end 
            loopcount = loopcount + 1;  % Increase loop count by one
            for cntr2=1:size(y,2)      % Add the identified loop to the set
                fbl(loopcount,cntr2) = y(cntr1,cntr2);
            end
            % Subtract the minimum elasticity from the elasticity
            % value of each of the transitions of the identified loop
            for cntr2=1:bg-1           
                E(y(cntr1,cntr2+1),y(cntr1,cntr2)) = E(y(cntr1,cntr2+1),y(cntr1,cntr2)) - min_elast;
            end
            E(y(cntr1,1),y(cntr1,bg)) = E(y(cntr1,1),y(cntr1,bg)) - min_elast;
            % Assign the minimum elasticity as the characteristic
            % elasticity of the identified loop
            els(loopcount) = min_elast;
        end % end for-cntr1
    end     % end for-i
end        % end for-bg

L.fbl=fbl; L.els=els;
  
    
function cycles = cycle(R)
%CYCLE Cycle partitions of reachability matrix.
%		CYCLES = cycle(R). The full list of elements in 
%       a cycle is reported in CYC.
%       Modified from Oliva's levels(R) function.

killc=sparse(size(R,1),1); 				% Constant to facilitate removing 
killr=sparse(1,size(R,2));  			% Rows and Columns

c = 1;
cycles=[];
control=sparse(ones(1,size(R,2)));		% All elements in the matrix
C=sparse(1,size(R,2));

while (nnz(control)>0)                  % While there are elements left
    T=sparse(1,size(R,2));
    for i=find(control)					% for all elements
        P = R(:,i)';					% id predecesors
        if (nnz(P)>0)					% if are any left
            S = R(i,:);					% id successors
            if (nnz(S&(~(P&S)))==0)		% if no e of S that are ~e of P&S
                T(1,i)=1;               %(i.e. if elements in P or S ath the same or lower level)
                if (nnz(P&S)>1)			% If at least one element both a predecessor and successor 
                    cycles(c,:)=S;         % We have a cycle: Place it in matrix
                    fprintf('Cycle is %d\n',c);
                    c=c+1;              % Increment cycle counter
                    S(1,i)=0;			% Eliminate original element from loop set
                    for m=find(S)		% Remove cycle from R matrix
                        R(m,:)=killr;
                        R(:,m)=killc;
                    end
                    C=C|S;
                end
            end
        end
    end
    for m=find(T)		% Remove elements sharing the same level from R matrix
        R(m,:)=killr;
        R(:,m)=killc;
    end
    control = control&(~(T|C));
end
cycles=sparse(cycles);

function y=loops(a,sz)
% LOOPS(a,sz)	Identifies directed cycles of size sz within a maximal
%               cycle set, i.e., a universal reachability matrix
%
%               Modified from Oliva's loops(R) function.
%               Pgs. 336-341 Societal Systems (Warfield 1989)

y=[];                                       %initialize output
D=vrtx_dist(a);
for i=1:5
% fprintf('Distance matrix %d %d %d %d %d %d\n', D(i,1), D(i,2), D(i,3), D(i,4), D(i,5), i)
end
GP=tril(D)'+triu(D);						% Length of cycle by adding paths a->b and b->a
% GP=GP - diag(diag(a));                      % Resize the length of self-circuits to 1
for i=1:5
% fprintf('Length matrix %d %d %d %d %d %d\n', GP(i,1), GP(i,2), GP(i,3), GP(i,4), GP(i,5), i)
end
maxloop=max(max(GP));
% fprintf('\n What is MAX loop? %d\n\n', maxloop);

k=1;

[r,c]=find(GP==1);                          % identify self loops
for i=1:size(r,1)                           % for each identified self loop 
    y(k,1)=r(i);                            % add it to the loop set
    k=k+1;                                  % increase counter by one
end

i = sz;                       			% for all potential loops whose size is equal to sz
[r,c]=find(GP==i);						% id elements involved in such loops
for j=1:size(r,1)   					% for all elements of loops same size
    if (r(j)~=0)    					% Check if loop has not been cleared
        loop=track_vrtx(r(j),c(j),D);
        for m=find(ismember(r(j+1:size(r,1)),loop))'	% for the roots of loops of same length in loop
            if size(m,1)>0								% extra if added for v 6.5 compatibility  -- in old ver. empty m would not enter for
                if ~isempty(find(loop==c(j+m)))				% if pair is already id
                    r(j+m)=0;								% remove from list
                end
            end
        end
        if size(loop,2) == size(unique(loop),2) 		% enter only if not figure 8
            y(k,1:size(loop,2))=loop;	
            k=k+1;
        end
    end
end

function y = vrtx_dist(a)
% DIST(m).		Calculates the reachability distance matrix
%				from a binary threshold matrix.
%
%               Function with the same name in MSA of Oliva.
%				Pg. 336 Societal Systems (Warfield 1989)

A=bsum(a,eye(size(a)));  	%Identity added to initial matrix
G = a;                      % Initialize the distance matrix
Bprev=A;
for i=2:size(a,1)			%Formula Pg. 336
	Bcur=bprd(A,Bprev);
	G=G+i*(Bcur-Bprev);
	Bprev=Bcur;
end
y=G;

function y=track_vrtx(u,v,dmatrix)
% The track-vrtx function derives "one" of the loops from vertex
% u through vertex v in a distance matrix.
% Function with the same name in MSA of Oliva.
% Pg. 340 Societal Systems (Warfield, 1989)

k=1;	x(k)=u;									% The u->v path
for j=1:(dmatrix(u,v)-1) 						% for 1 to (lenght of u->v path)-1
	rh=find(dmatrix(:,v) == (dmatrix(u,v)-j))'; % pred of v with step=length(u->v)-j
	lh=find(dmatrix(x(k),:) == 1); 				% succ of last element in path step=1			
	k=k+1;		
	for i=1:size(rh,2)
		if (find(lh == rh(1,i)))
			x(k)=rh(1,i);
			break;
        end
    end
end
								
k=k+1;	x(k)=v;									% The v->u path
for j=1:(dmatrix(v,u)-1)	 					% for 1 to (lenght of v->u path)-1
	rh=find(dmatrix(:,u) == (dmatrix(v,u)-j))';	% pred of u with step=length(v->u)-j
	lh=find(dmatrix(x(k),:) == 1); 				% succ of last element in path step=1			
	k=k+1;
	for i=1:size(lh,2)
		if (find(rh == lh(1,i)))
			x(k)=lh(1,i);
			break;
        end
    end
end
y=x;

function y=ils(f,a,order)
% ILS	Independent Loop Set.
%		ILS(f,a) reduces a set of feedback loops in a system to the maximal
%		number of independent loops. The f matrix has each loop in a row
%		-- assending in lenght -- with the elements listed across in columns.
%		The a matrix is the full adjaciancy matrix.
%       Function with the same name in MSA of Oliva.

b=(sparse(zeros(size(a))));				% matrix w/zeros to control base
k=1;									% counter for output matrix
L=ones(1,size(f,1));                    % intialize loop set to contain all loops

E(size(f,1))=struct('edges',sparse(zeros(size(a))));                 % initialize space for loop edges
for i=1:size(f,1)						% for each loop i in f
	length=nnz(f(i,:));						% identify its length
	T=sparse(zeros(size(a)));
	for j=1:length-1						% for each vertex in loop i
		T(f(i,j),f(i,j+1))=1;					% id edge between edge and next
    end
	T(f(i,length),f(i,1))=1;				% id edge from last to first node
	E(i).edges=T;
end
                                
while nnz(L)>0							% while loops in length set
	for j=find(L)						    % for every loop left
		c(j)=nnz((E(j).edges-b)>0); 			% id # of new edges contributed by loop
        if (c(j)==0)             % if loop does not contribute remove it from set
            L(j)=0;
        end			    
    end
	if order==-1 
        m=max(nonzeros(c));	    % if inverse order id maximal contribution           
    else 
        m=min(nonzeros(c)); 
    end			% else id minimal contribution
	if (m>0)                                % if at least one contribution
		m=find(c==m);                           % id loops with extreme contribution
        if order==-1  
            m=fliplr(m); 
        end         % flip so longest loop is first
		m=m(1);                                 % truncate pointer to the fist loop (shortest or longest)
		b=bsum(b,E(m).edges);						% add edges to base
		l(k,:)=f(m,:);							% add loop to output list
% 		fprintf(' k is: %d\n', k);
        k=k+1;									% increment loop counter
		L(m)=0;									% remove loop from length set
		c(m)=0;                                 % clear contribution vector
    end
end

[r,c]=find(l);		
l=l(1:max(r),1:max(c));					% clip the matrix to minimum size 
										% sort output by length
for i=1:size(l,1) 						% for all loops on list
	l_length(i)=nnz(l(i,:)); 				% id their length
end
k=1;									% initialize output counter
for i=sort(l_length)				% for ordered length
	j=find(l_length==i);					% find loop
	j=j(1);									% truncate
	y(k,:)=l(j,:);							% exchange
	k=k+1;
	l_length(j)=0;
end
function r=reach(i)
% REACH 	RCH=reach(ADJ)
%	    	Derives the full reachability matrix from an adjacency 
%			matrix (ADJ).
%			If ADJ is the structured array containing the adjacency 
%			matrix (.adj) and the vector with variable names (.nms), 
%			reach returns RCH with the same structure (.adj & .nsm),
%			but .adj containing the reachability matrix. 

% Adds the Identity matrix and raises to the max. distance possible
% The function uses the binary algebraic operations of bsum and bpwr.
% Function with the same name in MSA of Oliva.
% Pg. 233 Societal Systems (Warfield 1989)

if isa(i,'struct') 
	w=i.adj; 
else
	w=i;
end

if (size(w,1) ~= size(w,2)) 
    error('Argument is not a square matrix'); 
end
	
W=bsum(w,speye(size(w))); 
S=W;
P=bpwr(W,2);
while ((any(any(P>S))) ~= 0)
	S=P;
	P=bpwr(P,2);
end
if isa(i,'struct')
	r.adj=P; 
	r.nms=i.nms;
else
	r=P;
end

function y=bsum(a,b)
% BSUM Y=BSUM(A,B)
%	   Performs the binary addition of two matrices.
%	   Matrices have to be the same size.
%      Function with the same name in MSA of Oliva.

y=(a+b);
c=find(y>1);			%id nonzero elements that are not one	
y(c)=ones(size(c));		%truncate back to ones


function y=bprd(a,b)
% BPRD Y=PBRD(A,B)
%	   Performs the binary product of two matrices.
%	   The number of rows on the first factor has to match the 
%	   number of columns on the second.
%      Function with the same name in MSA of Oliva.

y =(a*b);
c=find(y>1);			%id nonzero elements that are not one	
y(c)=ones(size(c));		%truncate back to ones

function y=bpwr(a,exp)
% BPWR Y=BPRW(A,EXP) 
%	   Performs the power function in binary algebra (A^EXP).
%      Function with the same name in MSA of Oliva.

y=(a^exp);
c=find(y>1);	%id nonzero elements that are not one	
y(c)=ones(size(c));		%truncate back to ones
