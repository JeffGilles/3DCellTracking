function I = filter(this, varargin)
%Tracker.filter Filters trajectories
%

% === Input ===============================================================

p = inputParser;
addParameter(p, 'numel', [], @isnumeric);
parse(p, varargin{:});

fnumel = p.Results.numel;

% =========================================================================

% --- Number of elements

if ~isempty(fnumel)
    
    this.filters.numel = fnumel;
    
    n = arrayfun(@(x) numel(x.t), this.traj);
    Inumel = n>=fnumel(1) & n<=fnumel(2);
    
end

% --- Output --------------------------------------------------------------

I = find(Inumel);

if ~nargout
    
    this.traj = this.traj(I);
    
end