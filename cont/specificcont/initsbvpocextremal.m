function ocTrj=initsbvpocextremal(ocObj,solObj,varargin)
ocTrj=octrajectory();
opt=[];
if nargin>=4
    opt=varargin{1};
end
if isempty(opt)
    opt=defaultocoptions;
end
ZeroDeviationTolerance=getocoptions(opt,'GENERAL','ZeroDeviationTolerance');
TrivialArcMeshNum=getocoptions(opt,'GENERAL','TrivialArcMeshNum');

if isdynprimitive(solObj)
elseif isoctrajectory(solObj)
    ocTrj=octrajectory(solObj); % if solObj is an ocasymptotic the octrajectory part is returned
    vinfo=violationinfo(ocTrj);
    arcpos=arcposition(ocTrj);
    toalarcnum=arcnum(ocTrj);
    %inftimetrans=inftimetransformation(ocTrj);
    if ~isempty(vinfo)
        for ii=numel(vinfo):-1:1 % running backwards through the loop guarantees that the actual arc number, that has to be changed, is not influenced by previous violation adaptations
            newarcinfo(ii)=determinenewarc(ocObj,vinfo(ii));
            actarcnum=vinfo(ii).arcnum;
            actarcposition=arcpos(:,actarcnum);
            for jj=size(newarcinfo(ii).violationposition,2):-1:1
                if newarcinfo(ii).violationposition(1,jj)==actarcposition(1)
                    % violation occurs at the beginning of the arc
                    if ~(actarcnum==1)
                        % this case is an indication that the actual
                        % arc shrinks to a point, i.e., the arc can be
                        % removed
                        dintv=diff(arcinterval(ocEx));
                        if abs(dintv(actarcnum))<ZeroDeviationTolerance || dintv(actarcnum)<0 && actarcnum<toalarcnum
                            ocTrj=removearc(ocTrj,actarcnum,actarcposition);
                            toalarcnum=arcnum(ocTrj);
                        else
                            ocmatmsg('A violation at the beginning of an arc, different from the first/last arc, usually indicates that the entire arc can be romved./nThe actual result is different from that situation and can therefore not be processed automatically.\n')
                            return
                        end
                    elseif actarcnum==1
                        % violation occurs at the beginning of the first arc
                        ocTrj=addarc(ocTrj,actarcnum,actarcposition,vinfo(ii).constraintvalue(constraintindex(ii),:),newarcinfo(ii).violationposition(:,jj),TrivialArcMeshNum,'l');
                        toalarcnum=arcnum(ocTrj);
                    end
                elseif newarcinfo(ii).violationposition(2,jj)==actarcposition(2)
                    % violation ends at the end of the arc
                    if isequilibrium(solObj)
                        ocmatmsg('Violation ends at the end of the arc.\nSince the ocasymptotic converges to an equilibrium this result refers to a numerical artefact.\Assure that the path ends near the equilibrium, maybe increase time horizon.\n')
                        return
                    end
                else
                    % violation occurs inbetween the arc
                    ocTrj=addarc(ocTrj,actarcnum,vinfo(ii).constraintvalue(newarcinfo(ii).constraintindex(jj),:),newarcinfo(ii).violationposition(:,jj),newarcinfo(ii).newarcarg(jj),TrivialArcMeshNum,'m');
                end
            end
        end
    end
end

ocTrj=postprocess(ocTrj); % increase number of grid points for each arc to at least 10

if isocasymptotic(solObj)
    ocTrj=ocasymptotic(ocTrj,limitset(solObj));
end

function ocTrj=removearc(ocTrj,actarcnum,actarcposition)

toalarcnum=arcnum(ocTrj);
arcpos=arcposition(ocTrj);
removepos=actarcposition(1):actarcposition(2);

% reduce the time interval right to the removed arc by one
ocTrj.x(actarcposition(2)+1:arcpos(2,toalarcnum))=ocTrj.x(actarcposition(2)+1:arcpos(2,toalarcnum))-1;

% remove the arc
ocTrj.arcinterval(actarcnum+1)=[];
ocTrj.x(removepos)=[];


function ocTrj=addarc(ocTrj,actarcnum,constraintvalue,violationposition,newarcarg,N,addflag)
toalarcnum=arcnum(ocTrj);
arcpos=arcposition(ocTrj);
actarcpos=arcpos(1,actarcnum):arcpos(2,actarcnum);
actviolationpos=violationposition(1):violationposition(2);
numviolationpt=numel(actviolationpos);
t=ocTrj.x(actarcpos);
y=ocTrj.y(:,actarcpos);
switch addflag
    case 'l'
        % add arc starting at the left side of the interval
        if numviolationpt==1
            ynew=[repmat(y(:,actviolationpos),1,N) y];
            tnew=[repmat(t(actviolationpos),1,N) t+1];
        else
        end
    case 'r'
        % add arc ending at the right side of the interval
    case 'm'
        % add arc lying inbetween the interval
        if numviolationpt==1
        else
            % extended time interval including the grid point lying on the
            % left and right of the interval where the violation occurs
            txtended=t([violationposition(1)-1 actviolationpos violationposition(2)+1]);
            yxtended=y(:,[violationposition(1)-1 actviolationpos violationposition(2)+1]);
            cxtended=constraintvalue([violationposition(1)-1 actviolationpos violationposition(2)+1]);
            % determine intersection points on the time interval where
            % constraint value becomes zero
            cl=[cxtended(1) cxtended(2)];
            cr=[cxtended(numviolationpt+1) cxtended(numviolationpt+2)];
            tl=[txtended(1) txtended(2)];
            tr=[txtended(numviolationpt+1) txtended(numviolationpt+2)];
            tzerol=(tl(1)*cl(2)-tl(2)*cl(1))/(cl(2)-cl(1));
            tzeror=(tr(1)*cr(2)-tr(2)*cr(1))/(cr(2)-cr(1));
            tint=linspace(tzerol,tzeror,N);
            yint=interp1(txtended,yxtended.',tint).';
            s1=[t(1:violationposition(1)-1) tzerol]/(tzerol-t(1))+actarcnum-1;
            s2=(tint-tint(1))/(tint(N)-tint(1))+actarcnum;
            s3=([tzeror t(violationposition(2)+1:end)]-tzeror)/(t(end)-tzeror)+actarcnum+1;
            y1=[y(:,1:violationposition(1)-1) yint(:,1)];
            y2=yint;
            y3=[yint(:,end) y(:,violationposition(2)+1:end)];
            ocTrj.arcinterval=[ocTrj.arcinterval(1:actarcnum) tzerol tzeror ocTrj.arcinterval(actarcnum+1)];
            ocTrj.arcarg=[ocTrj.arcarg(1:actarcnum) newarcarg ocTrj.arcarg(actarcnum)];
            ocTrj.x=[ocTrj.x(1:arcpos(1,actarcnum)-1) s1 s2 s3 ocTrj.x(arcpos(2,actarcnum)+1:arcpos(2,toalarcnum))];
            ocTrj.y=[ocTrj.y(:,1:arcpos(1,actarcnum)-1) y1 y2 y3 ocTrj.y(:,arcpos(2,actarcnum)+1:arcpos(2,toalarcnum))];
            idx=find(diff(ocTrj.x)==0);
            ocTrj.arcposition=[[1 idx+1];[idx numel(ocTrj.x)]];
        end
end

function ocTrj=postprocess(ocTrj)
n=10;

% remove solver information that is not actual any more, since the grid and data has changed
datadependentfields={'coeff','tangent','continuationparameter'};
ocTrj.solver='';
solverinfoStruct=solverinfo(ocTrj);
fn=fieldnames(solverinfoStruct);
for ii=1:numel(fn)
    if any(strcmp(datadependentfields,fn{ii}))
        solverinfoStruct=rmfield(solverinfoStruct,fn{ii});
    end
end
ocTrj.solverinfo=solverinfoStruct;

fidx=find(diff(ocTrj.arcposition)<n);

if isempty(fidx)
    return
end
narc=arcnum(ocTrj);
nfidx=numel(fidx);
if fidx(nfidx)<narc
    tnew=[ocTrj.x(ocTrj.arcposition(1,fidx(nfidx)+1):ocTrj.arcposition(2,narc))];
    ynew=[ocTrj.y(:,ocTrj.arcposition(1,fidx(nfidx)+1):ocTrj.arcposition(2,narc))];
else
    tnew=[];
    ynew=[];
end
for ii=nfidx:-1:1
    actposition=ocTrj.arcposition(:,fidx(ii));
    tint=linspace(ocTrj.x(actposition(1)),ocTrj.x(actposition(2)),n);
    yint=interp1(ocTrj.x(actposition(1):actposition(2)),ocTrj.y(:,actposition(1):actposition(2)).',tint).';
    tnew=[tint tnew];
    ynew=[yint ynew];
end
ocTrj.x=tnew;
ocTrj.y=ynew;
idx=find(diff(ocTrj.x)==0);
ocTrj.arcposition=[[1 idx+1];[idx numel(ocTrj.x)]];


