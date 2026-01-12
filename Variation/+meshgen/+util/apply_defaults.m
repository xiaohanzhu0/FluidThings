function s = apply_defaults(s, defaults)
%APPLY_DEFAULTS Fill missing or empty fields in a struct.
    if nargin < 1 || isempty(s)
        s = struct();
    end
    if nargin < 2 || isempty(defaults)
        return;
    end
    fields = fieldnames(defaults);
    for i = 1:numel(fields)
        name = fields{i};
        if ~isfield(s, name) || isempty(s.(name))
            s.(name) = defaults.(name);
        end
    end
end
