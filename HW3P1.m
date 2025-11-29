%cvx_solver sdpt3 %default, used in project
%%
clear;
clc;
fileID = fopen('TSP_1000_euclidianDistance.txt','r'); %"euclidian" or "random"
len = fscanf(fileID, '%d');
fscanf(fileID, '%s %s %s', 3);
formatSpec = '%d %d %f';
sizeA = [3 Inf];
ABC = transpose(fscanf(fileID,formatSpec,sizeA));

A = ABC(:, 1);
B = ABC(:, 2);

C = ones(len, len) * 65535;

for i = 1:size(A, 1)
    C(A(i), B(i)) = ABC(i, 3);
end

upper = triu(ones(len, len), 1);
inGraph = ones(len, 1);

%%
cvx_begin
    cvx_solver gurobi
    params.TimeLimit = 60;
    
    variable x(len, len) binary;
    variable y(len, len) nonnegative;

    minimize sum(sum(x .* C))
        subject to
            x <= upper;
            sum(x + transpose(x), 2) == 2; 
            ((len - 1)/2)*x(1, 2:end) + sum(y(:, 2:end)) == sum(transpose(y(2:end, :))) + 1;
            (y(2:end, 2:end) + transpose(y(2:end, 2:end)) .* upper(2:end, 2:end)) == (((len - 3)/2)*x(2:end, 2:end));
cvx_end

xM = full(x);
G = graph(xM + transpose(xM));
D = dfsearch(G, 1);
DS = size(D, 1);
loops = 0;

while (DS < len)
    loops = loops + 1;
    D1max = 0;
    D2max = 0;
    D11 = D(1);
    D12 = D(end);
    
    for i = 1:1:(DS - 1)
        ind = D(i);
        inGraph(ind) = 0;
        Dtemp = D(i + 1);
        
        if C(min(ind,Dtemp), max(ind, Dtemp)) > D1max;
            D11 = min(ind, Dtemp);
            D12 = max(ind, Dtemp);
            D1max = C(min(ind,Dtemp), max(ind, Dtemp));
        end
    end

    ind = D(end);
    inGraph(ind) = 0;
    Dtemp = D(1);
    
    if C(Dtemp, ind) > D1max;
        D11 = Dtemp;
        D12 = ind;
        D1max = C(Dtemp, ind);
    end

    ind = find(inGraph == 1, 1, 'first');
    D2 = dfsearch(G, ind);
    D21 = D2(1);
    D22 = D2(end);
    
    for i = 1:1:(size(D2, 1) - 1)
        ind = D2(i);
        Dtemp = D2(i + 1);
        
        if C(min(ind,Dtemp), max(ind, Dtemp)) > D2max;
            D21 = min(ind, Dtemp);
            D22 = max(ind, Dtemp);
            D2max = C(min(ind,Dtemp), max(ind, Dtemp));
        end
    end

    ind = D2(end);
    Dtemp = D2(1);
    
    if C(Dtemp, ind) > D2max;
        D21 = Dtemp;
        D22 = ind;
        D2max = C(Dtemp, ind);
    end

    xM(D11, D12) = 0;
    xM(D21, D22) = 0;
    

    if (C(min(D11, D22), max(D11, D22)) + C(min(D21, D12), max(D21, D12))) < (C(min(D11, D21), max(D11, D21)) + C(min(D22, D12), max(D22, D12))) 
        xM(min(D11, D21), max(D11, D21)) = 1;
        xM(min(D22, D12), max(D22, D12)) = 1;
    else
        xM(min(D11, D22), max(D11, D22)) = 1;
        xM(min(D21, D12), max(D21, D12)) = 1;
    end

    G = graph(xM + transpose(xM));
    D = dfsearch(G, 1);
    DS = size(D, 1);
end

fprintf("\n\n");

for ind = 1:1:len-1
    fprintf("%d, ", D(ind));
end

fprintf("%d\n", D(len))

final = sum(sum(xM .* C));