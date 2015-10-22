function [K0time K1time K2time K3time] = PageRankPipeline(SCALE,EdgesPerVertex,Nfile);

  Nmax = 2.^SCALE;                                 % Max vertex ID.
  M = EdgesPerVertex .* Nmax;                      % Total number of edges.
  myFiles = 1:Nfile;                               % Set list of files.
  %myFiles = global_ind(zeros(Nfile,1,map([Np 1],{},0:Np-1)));   % PARALLEL.
  tab = char(9);
  nl = char(10);
  Niter = 20;                                      % Number of PageRank iterations.
  c = 0.15;                                        % PageRank damping factor.

  disp(['Number of Edges: ' num2str(M) ', Maximum Possible Vertex: ' num2str(Nmax)]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kernel 0: Generate a Graph500 Kronecker graph and save to data files.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Kernel 0: Generate Graph, Write Edges']);
  tic;

    for i = myFiles
      fname = ['data/K0/' num2str(i) '.tsv'];
      disp(['  Writing: ' fname]);  % Read filename.
      rand('seed',i);                              % Set random seed to be unique for this file.
      [u v] = KronGraph500NoPerm(SCALE,EdgesPerVertex./Nfile);      % Generate data.
 
      edgeStr = sprintf(['%16.16g' tab '%16.16g' nl],[u'; v']);    % Convert edges to strings.
      edgeStr = edgeStr(edgeStr ~= ' ');                            % Remove extra spaces.
      StrFileWrite(edgeStr,fname);                                  % Write string to file.
    end

  K0time = toc;
  disp(['K0 Time: ' num2str(K0time) ', Edges/sec: ' num2str(M./K0time)]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kernel 1: Read data, sort data, and save to files.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Kernel 1: Read, Sort, Write Edges']);
  edgeStr = '';
  tic;

    % Read in all the files into one array.
    for i = myFiles
      fname = ['data/K0/' num2str(i) '.tsv'];
      disp(['  Reading: ' fname]);  % Read filename.
      edgeStr = [edgeStr StrFileRead(fname)];
    end

    % Sort by start edge.
    uv = sscanf(edgeStr,'%f');                      % Convert string to numeric data.
    u = uv(1:2:end);                                % Get starting vertices.
    v = uv(2:2:end);                                % Get ending vertices.
    [u sortIndex] = sort(u);                        % Sort starting vertices.
    v = v(sortIndex);                               % Carry ending vertices.

    % Write all the data to files.
    j = 1;                                                             % Initialize file counter.
    for i = myFiles
      jEdgeStart = ((j-1).*(size(u,1)./numel(myFiles))+1);            % Compute first edge of file.
      jEdgeEnd = ((j).*(size(u,1)./numel(myFiles)));                  % Compute last edge of file.
      uu = u(jEdgeStart:jEdgeEnd);                                   % Select start vertices.
      vv = v(jEdgeStart:jEdgeEnd);                                   % Select end vertices.
      fname = ['data/K1/' num2str(i) '.tsv'];
      disp(['  Writing: ' fname]);  % Create filename.
      edgeStr = sprintf(['%16.16g' tab '%16.16g' nl],[uu'; vv']);   % Convert edges to strings.
      edgeStr = edgeStr(edgeStr ~= ' ');                            % Remove extra spaces.
      StrFileWrite(edgeStr,fname);                                   % Write string to file.
      j = j + 1;                                                       % Increment file counter.
    end

  K1time = toc;
  disp(['K1 Time: ' num2str(K1time) ', Edges/sec: ' num2str(M./K1time)]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kernel 2: Read data, filter data.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Kernel 2: Read, Filter Edges']);
  edgeStr = '';
  tic;

    % Read in all the files into one array.
    for i = myFiles
      fname = ['data/K1/' num2str(i) '.tsv'];
      disp(['  Reading: ' fname]);  % Read filename.
      edgeStr = [edgeStr StrFileRead(fname)];
    end

    % Construct adjacency matrix.
    uv = sscanf(edgeStr,'%f');       % Convert string to numeric data.
    u = uv(1:2:end);                 % Get starting vertices.
    v = uv(2:2:end);                 % Get ending vertices.
    A = sparse(u,v,1,Nmax,Nmax);     % Create adjacency matrix.

    % Filter and weight the adjacency matrix.
    din = sum(A,1);                    % Compute in degree.
    A(:,din == max(din)) = 0;          % Eliminate the super-node.
    A(:,din == 1) = 0;                 % Eliminate the leaf nodes.
    dout = sum(A,2);                   % Compute the out degree.
    i = dout > 0;                      % Find vertices with outgoing edges.
    DoutInv = sparse(find(i),find(i),1./dout(i),Nmax,Nmax);   % Create diagonal weight matrix.
    A = DoutInv * A;                   % Apply weight matrix.

  K2time = toc;
  disp(['K2 Time: ' num2str(K2time) ', Edges/sec: ' num2str(M./K2time)]);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Kernel 3: Compute PageRank.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['Kernel 3: PageRank']);
  tic;

    r = rand(Nmax,1);                     % Generate a random starting rank.
    r = r ./ norm(r,1);                   % Normalize
    a = ones(Nmax,1) .* (1-c) ./ Nmax;    % Create damping vector

    for i=1:Niter
      r = A * (r .* c) + a;               % Compute PageRank.
    end

  K3time = toc;
  disp(['  Sum of PageRank: ' num2str(sum(r))]);     % Force all computations to occur.
  disp(['K3 Time: ' num2str(K3time) ', Edges/sec: ' num2str(Niter.*M./K3time)]);

%whos
%keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PageRank Pipeline Benchmark
% Architect: Dr. Jeremy Kepner (kepner@ll.mit.edu)
% MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) <2015> Massachusetts Institute of Technology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

