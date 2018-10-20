function HHMMSS=secs2str(t)
% It converts numerical seconds in printable HH:MM:SS format
if ~nargin
    help secs2str
    return
end

HH=fix(t/(60*60)); t=t-HH*60*60; HH=num2str(HH); HH=[num2str(zeros(1,2-numel(HH))) HH];
MM=fix(t/(60));    t=t-MM*60;    MM=num2str(MM); MM=[num2str(zeros(1,2-numel(MM))) MM];
SS=fix(t);                       SS=num2str(SS); SS=[num2str(zeros(1,2-numel(SS))) SS];

% HHMMSS=[HH ':' MM  ':' SS];
HHMMSS=[HH 'hr ' MM  'min ' SS 'sec'];
