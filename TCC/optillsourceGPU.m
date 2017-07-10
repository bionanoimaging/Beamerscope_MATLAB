function [ x_opt, fval, history ] = optillsourceGPU( method, eigenvalue, eigenfunction, OPTparams, object )
%optillsource Summary of this function goes here
%   Detailed explanation goes here

geigenvalue = gpuArray(eigenvalue);
geigenfunction = gpuArray(eigenfunction);
gobject = gpuArray(object);

% Set up shared variables with OUTFUN
history.x = [];
history.fval = [];
history.meshsize = [];

if( strcmp(OPTparams.optmethod, 'swarm'))
    
    if isfield(OPTparams, 'xmask')
        xmask = OPTparams.xmask;
    else
        % if no mask available use default
        xmask = ones(1,36);
    end
    
    if isfield(OPTparams, 'x0')
        x0 = OPTparams.x0(find(xmask));
    else
        % reduce number of used kernels to 2
        x0 = zeros(1,36);
    end
    
    
    
    
    nvars = size(x0,2);
    
    
    if isfield(OPTparams, 'lb')
        lb = OPTparams.lb;
        lb = lb*ones(1,nvars);
    else
        lb = zeros(1,nvars);
        if strcmp(method, 'Zernike')
            lb = ones(1,nvars)*(-10);
        end
    end
    
    if isfield(OPTparams, 'ub')
        ub = OPTparams.ub;
        ub = ub*ones(1,nvars);
    else
        ub = ones(1,nvars);
        if strcmp(method, 'Zernike')
            ub = ones(1,nvars)*(10);
        end
    end
    
    if isfield(OPTparams, 'A')
        A = OPTparams.A;
    else
        A = [];
    end
    
    if isfield(OPTparams, 'b')
        b = OPTparams.b;
    else
        b = [];
    end
    
    if isfield(OPTparams, 'Aeq')
        Aeq = OPTparams.Aeq;
    else
        Aeq = [];
    end
    
    if isfield(OPTparams, 'beq')
        beq = OPTparams.beq;
    else
        beq = [];
    end
    
    if isfield(OPTparams, 'MaxIter')
        MaxIter = OPTparams.MaxIter;
    else
        MaxIter = 10;
    end
    
    if isfield(OPTparams, 'kernelOrder')
        kernelOrder = OPTparams.kernelOrder;
    else
        % reduce number of used kernels to 2
        kernelOrder = 2;
    end
    
    if isfield(OPTparams, 'MaxTime')
        MaxTime = OPTparams.MaxTime;
    else
        % reduce number of used kernels to 2
        MaxTime = 60;
    end
    
    if isfield(OPTparams, 'StallIterLimit')
        StallIterLimit = OPTparams.StallIterLimit;
    else
        % reduce number of used kernels to 2
        StallIterLimit = 60;
    end
    
    
    
    % obtain spectrum of object-transmission function
    objectspectrum=fftshift(ifft2(ifftshift(object)));
    
    % assign objective function
    fun = @psobjTCC;
    objective = @(x)fun(x, xmask, objectspectrum, geigenfunction, geigenvalue, kernelOrder, object, OPTparams.qualitymethod);
    
    % initiliaze optimization parameters
    %, 'MaxTime', 60
    
    % swarmsize suggestion: https://www.researchgate.net/post/On_what_basis_do_we_select_the_swarm_size_for_any_application_in_Particle_Swarm_Optimization_Does_it_vary_with_the_type_of_PSO_used
    options=optimoptions('particleswarm','SwarmSize',nvars*15,'Display','iter','UseParallel',true,  'MaxTime', MaxTime, 'OutputFcn',@outfunNest); %'PlotFcns',@pswplotbestf
    options.InitialSwarm = (x0);
    options.StallIterLimit = StallIterLimit;
    options.MaxIter =  MaxIter;
    options.TolFun = 1e-8;
    
    % do optimization
    [x_opt_temp,fval,exitflag,output_opt] = particleswarm(objective, nvars, lb, ub, options);%,A,b,Aeq,beq,lb,ub,options)
    
    % fit the variables back in the x-vector/coefficiencts vector
    x_opt = zeros(size(OPTparams.x0));
    x_opt(find(xmask)) = x_opt_temp;
    
    for(i = 1:size(history.x,1))
        
        history_out = zeros(size(OPTparams.x0));
        history_out(find(xmask)) = history.x(i,:);
        histroy_temp(:,i) = history_out;
    end
    
    history.x = histroy_temp;
    
    
elseif( strcmp(OPTparams.optmethod, 'grid'))
    
    
    % setup the optimization parameters
    
    
    
    
    if isfield(OPTparams, 'xmask')
        xmask = OPTparams.xmask;
    else
        % if no mask available use default
        xmask = ones(1,36);
    end
    
    if isfield(OPTparams, 'x0')
        find(xmask)
        x0 = OPTparams.x0(find(xmask));
    else
        % reduce number of used kernels to 2
        x0 = zeros(1,36);
    end
    
    
    
    
    nvars = size(x0,2);
    
    
    if isfield(OPTparams, 'lb')
        lb = OPTparams.lb;
        lb = lb*ones(1,nvars);
    else
        lb = zeros(1,nvars);
        if strcmp(method, 'Zernike')
            lb = ones(1,nvars)*(-10);
        end
    end
    
    if isfield(OPTparams, 'ub')
        ub = OPTparams.ub;
        ub = ub*ones(1,nvars);
    else
        ub = ones(1,nvars);
        if strcmp(method, 'Zernike')
            ub = ones(1,nvars)*(10);
        end
    end
    
    if isfield(OPTparams, 'A')
        A = OPTparams.A;
    else
        A = [];
    end
    
    if isfield(OPTparams, 'b')
        b = OPTparams.b;
    else
        b = [];
    end
    
    if isfield(OPTparams, 'Aeq')
        Aeq = OPTparams.Aeq;
    else
        Aeq = [];
    end
    
    if isfield(OPTparams, 'beq')
        beq = OPTparams.beq;
    else
        beq = [];
    end
    
    if isfield(OPTparams, 'MaxIter')
        MaxIter = OPTparams.MaxIter;
    else
        MaxIter = 10;
    end
    
    if isfield(OPTparams, 'kernelOrder')
        kernelOrder = OPTparams.kernelOrder;
    else
        % reduce number of used kernels to 2
        kernelOrder = 2;
    end
    
    
    
    % obtain spectrum of object-transmission function
    gobjectspectrum=fftshift(ifft2(ifftshift(gobject)));
    
    % assign objective function
    fun = @psobjTCC;
    objective = @(x)fun(x, xmask, gobjectspectrum, eigenfunction, eigenvalue, kernelOrder, object, OPTparams.qualitymethod);
        
    options = psoptimset('Display','diagnose','PlotFcns',@psplotbestf,'MaxIter', MaxIter, 'TimeLimit', 60, 'UseParallel', 'always', 'OutputFcn', @outfunNestGrid);
    [x_opt_temp,fval,exitflag,output] = patternsearch(objective,x0,A,b,Aeq,beq,lb,ub,options);
    
    
    % fit the variables back in the x-vector/coefficiencts vector
    x_opt = zeros(size(OPTparams.x0));
    x_opt(find(xmask)) = x_opt_temp;
    
    for(i = 1:size(history.x,1))
        
        history_out = zeros(size(OPTparams.x0));
        history_out(find(xmask)) = history.x(i,:);
        histroy_temp(:,i) = history_out;
    end
    
    history.x = histroy_temp;
    
else
    error('no optimization method has been chosen!')
    
    
end



    function stop = outfunNest(optimValues,state)
        stop = false;
        x= optimValues.bestx;
        val = optimValues.bestfval;
        
        switch state
            case 'init'
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; val];
                history.x = [history.x; x];
                % Concatenate current search direction with
                % searchdir.
            case 'done'
            otherwise
        end
    end

    function [stop,options,optchanged] = outfunNestGrid(optimvalues,options,flag,p1,p2)
        
        stop = false;
        x= optimvalues.x;
        val = optimvalues.fval;
        meshsize = optimvalues.meshsize;
        
        optchanged = false;
        
        
        switch flag
            case 'init'
            case 'iter'
                % Concatenate current point and objective function
                % value with history. x must be a row vector.
                history.fval = [history.fval; val];
                history.x = [history.x; x];
                history.meshsize = [history.meshsize; meshsize];
                % Concatenate current search direction with
                % searchdir.
            case 'done'
            otherwise
        end
    end


end

