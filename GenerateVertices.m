%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenerateVertices
% Owner: Jay McMahon (The University of Colorado Boulder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Description:
%
%  This function loads the shape model and assigns parameters in desired
%  format for other functions.
%
%%% Inputs:
%
%    - ShapeModelFile  : Shape model OBJ file name
%
%%% Outputs:
%
%    - Polygon         : Shape model information
%
%%% Assumptions/References:
%
%   - None
%
%%% Note:
%
%   - None
%
%%% Dependencies:
%
%   - None
%
%%% Call
%
%   - read_obj.m
%
%%% Called by
%
%   - Run_SlopeMap.m
%
%%% Modification History:
%
%   15Oct14   Jay McMahon    Delivery 1 to SPOC
%   24Feb12   Yu Takahashi   original version
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Polygon = GenerateVertices(ShapeModelFile)

% Load data from OBJ file
[FacetIndex,verts] = read_obj(ShapeModelFile);

num_vertex = length(verts(:,1));     % [n.d.] Number of the vertices
Vertex = [(1:num_vertex)' verts];    % [km]   Positions of the vertices

%%%%%%%%%%%%%%%%%%%%%%%%
%% -- Facets (All) -- %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Cross-product of the second vertex minus the first vertex and the third vertex minus the first vertex must points out of
% the facet

num_facet                   = length(FacetIndex(:,1)); % [n.d.]    Number of the facets
Facet                       = zeros(3,3,num_facet);    % [km]      Positions of the three vertices to make facets
FacetNormal                 = zeros(num_facet,3);      % [n.d.]    Surface normal
FacetCenter                 = zeros(num_facet,3);      % [km]      Location of the facet center
FacetVertex_1               = zeros(num_facet, 3);     % [km]      Positions of the first verteces
FacetVertex_2               = zeros(num_facet, 3);     % [km]      Positions of the second verteces
FacetVertex_3               = zeros(num_facet, 3);     % [km]      Positions of the third verteces

for ii = 1:num_facet
    
    for ff = 1:3
        
        Facet(ff,:,ii) = Vertex(FacetIndex(ii,ff),2:4);  % [km] Position of one of the vertices to make a facet
        
        if ff == 1
            
            FacetVertex_1(ii,:) = Vertex(FacetIndex(ii,ff),2:4);
            
        elseif ff == 2
            
            FacetVertex_2(ii,:) = Vertex(FacetIndex(ii,ff),2:4);
            
        elseif ff == 3
            
            FacetVertex_3(ii,:) = Vertex(FacetIndex(ii,ff),2:4);
            
        end % For if
        
    end % For ff
    
    D                 = FacetVertex_1(ii,:);
    E                 = FacetVertex_2(ii,:);
    F                 = FacetVertex_3(ii,:);
    FacetCenter(ii,:) = (D + E + F)/3;
    
    G          = E - D;
    H          = F - D;
    N          = cross(G, H);
    
    FacetNormal(ii,:) = N/norm(N);
    
end % For ii

%%%%%%%%%%%%%%%%%%%%%%%
%% -- Edges (All) -- %%
%%%%%%%%%%%%%%%%%%%%%%%

EdgesIndex_All = zeros(3*num_facet,2); % [n.d.] Three edges for each facet, and two vertices make an edge

for ii = 1:num_facet
    
   EdgesIndex_1 = [min([FacetIndex(ii,1), FacetIndex(ii,2)]), max([FacetIndex(ii,1), FacetIndex(ii,2)])];
   EdgesIndex_2 = [min([FacetIndex(ii,2), FacetIndex(ii,3)]), max([FacetIndex(ii,2), FacetIndex(ii,3)])];
   EdgesIndex_3 = [min([FacetIndex(ii,3), FacetIndex(ii,1)]), max([FacetIndex(ii,3), FacetIndex(ii,1)])];
   
   EdgesIndex_All(3*(ii-1)+1:3*(ii-1)+3,:)   =  [EdgesIndex_1; EdgesIndex_2; EdgesIndex_3]; % An edge is defined [smaller index number , larger index number]
   
end % For ii

num_edges_all  = length(EdgesIndex_All(:,1)); % [n.d.] Number of the edges (the same edge is counted twice)
EdgesIndex     = zeros(num_edges_all/2,2);    % [n.d.] Redundant edges exist, and two vertices make an edge
count          = 0;

for ii = 1:num_edges_all
    
    [row, ~] = find(EdgesIndex_All(:,1) == EdgesIndex_All(ii,1) & EdgesIndex_All(:,2) == EdgesIndex_All(ii,2));
    
    if ( ~isempty(row) ) && ( EdgesIndex_All(ii,1) ~= 0 ) % zero doesn't count
        
        count                    = count + 1;
        EdgesIndex(count,:)      = EdgesIndex_All(row(1),:);
        EdgesIndex_All(row(2),:) = [0, 0];
        
    end % For if
    
end % For ii

num_edges                  = length(EdgesIndex(:,1)); % [n.d.] Number of the edges
Edges                      = zeros(2,3,num_edges);    % [n.d.] Two vertices make an edge, and there is "num_edges" number of edges
Edges_1                    = zeros(num_edges,3);      % [km]   Positions of the first vertices
Edges_2                    = zeros(num_edges,3);      % [km]   Positions of the second vertices
EdgeFacetNumberDirection_1 = zeros(num_edges,2);
EdgeFacetNumberDirection_2 = zeros(num_edges,2);

for ii = 1:num_edges
   
    Edges(1,:,ii) = Vertex(EdgesIndex(ii,1),2:4); % [km] First row is the position of the vertex with smaller index value
    Edges(2,:,ii) = Vertex(EdgesIndex(ii,2),2:4); % [km] First row is the position of the vertex with larger index value
    
    Edges_1(ii,:) = Vertex(EdgesIndex(ii,1),2:4);
    Edges_2(ii,:) = Vertex(EdgesIndex(ii,2),2:4);
    
    IndexVertex_1 = EdgesIndex(ii,1);
    IndexVertex_2 = EdgesIndex(ii,2);
    
    count = 0;
    
    for ff = 1:num_facet
        
        edge_on = 0;
        
        if ( FacetIndex(ff,1) == IndexVertex_1 ) && ( FacetIndex(ff,2) == IndexVertex_2 )
        
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = 1;
            
        elseif ( FacetIndex(ff,2) == IndexVertex_1 ) && ( FacetIndex(ff,1) == IndexVertex_2 )
            
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = -1;
            
        elseif ( FacetIndex(ff,1) == IndexVertex_1 ) && ( FacetIndex(ff,3) == IndexVertex_2 )
         
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = -1;
            
        elseif ( FacetIndex(ff,3) == IndexVertex_1 ) && ( FacetIndex(ff,1) == IndexVertex_2 )
            
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = 1;
            
        elseif ( FacetIndex(ff,2) == IndexVertex_1 ) && ( FacetIndex(ff,3) == IndexVertex_2 )
            
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = 1;
            
        elseif ( FacetIndex(ff,3) == IndexVertex_1 ) && ( FacetIndex(ff,2) == IndexVertex_2 )
            
            edge_on       = 1;
            count         = count + 1;
            EdgeDirection = -1;
            
        end
        
        if edge_on == 1
            
            if count == 1
                
                EdgeFacetNumberDirection_1(ii,:) = [ff,EdgeDirection];
                
            elseif count == 2
                
                EdgeFacetNumberDirection_2(ii,:) = [ff,EdgeDirection];
                break
                
            end % For if
                        
        end % For if
            
    end % For if
    
end % For ii



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -- Assigne all the outputs to Polygon -- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Polygon.Vertex                      = Vertex;
Polygon.Edges                       = Edges;
Polygon.Edges_1                     = Edges_1;
Polygon.Edges_2                     = Edges_2;
Polygon.EdgesIndex                  = EdgesIndex;
Polygon.EdgeFacetNumberDirection_1  = EdgeFacetNumberDirection_1;
Polygon.EdgeFacetNumberDirection_2  = EdgeFacetNumberDirection_2;
Polygon.FacetIndex                  = FacetIndex;
Polygon.Facet                       = Facet;
Polygon.FacetVertex_1               = FacetVertex_1;
Polygon.FacetVertex_2               = FacetVertex_2;
Polygon.FacetVertex_3               = FacetVertex_3;
Polygon.FacetNormal                 = FacetNormal;
Polygon.FacetCenter                 = FacetCenter;
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%