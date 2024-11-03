#include "Instance.hpp"


void node::print(){
    string T="N_";
    if(isInitialDepot()) T="O_";
    else if(isFinalDepot()) T="D_";
    else T="N_";
    cout<<T<<_id<<"("<<_x<<", "<<_y<<", "<<")";
    
}

Instance::Instance(string filename)
{
    _name=filename;
    
    ifstream input(filename);
    if(!input.is_open()){
        throw runtime_error("Could not open file "+filename);
    }
    assert(input.is_open()); 

    string name = getData(input, "NAME"); 
    string type = getData(input, "TYPE"); 
    string comment = getData(input, "COMMENT"); 
    string dimension = getData(input, "DIMENSION");
    string edge_weight_type = getData(input, "EDGE_WEIGHT_TYPE");
    string edge_weight_format = getData(input, "EDGE_WEIGHT_FORMAT");

    cout << "Name: " << name << endl; 
    cout << "Type: " << type << endl;
    cout << "Comment: " << comment << endl;
    cout << "Dimension: " << dimension << endl;
    cout << "Edge Weight Type: " << edge_weight_type << endl;
    cout << "Edge Weight Format: " << edge_weight_format << endl;

    if(name=="n/a" || type=="n/a" || comment=="n/a" || dimension=="n/a" || edge_weight_type=="n/a"){
        throw runtime_error("Could not find the basic information in file "+filename);
    }

 
    _name=name;
    _nnodes=stoi(dimension);
    _Nodes=new node*[_nnodes];

    bool hasCoords=true;
    if(edge_weight_type=="EXPLICIT"){
        if(!gotoSection(input, "DISPLAY_DATA_SECTION")){
            hasCoords=false;
        }
    }
    else if(edge_weight_type=="EUC_2D" || edge_weight_type=="CEIL_2D" || edge_weight_type=="ATT" || edge_weight_type=="GEO"){
        if(!gotoSection(input, "NODE_COORD_SECTION")){
            throw runtime_error("Could not find NODE_COORD_SECTION in file "+filename);
        }
    }
    else{
        throw runtime_error("Edge weight type "+edge_weight_type+" is not implemented");
    }

    if(!hasCoords && edge_weight_type!="EXPLICIT"){
        throw runtime_error("DISPLAY_DATA_SECTION is not implemented");
    }

    if(hasCoords){
        for(int i=0;i<_nnodes;i++){
            int id;
            double x,y;
            input>>id>>x>>y;
            _Nodes[i]=new node(id,x,y);
        }
    }
    else{
        for(int i=0;i<_nnodes;i++)
            _Nodes[i]=new node(i+1);
    }

    _D=NULL;
    if(edge_weight_type=="EXPLICIT"){
        readExplicit(input, edge_weight_format);
        _DM=(DistanceMode)DistanceMode::Explicit;
    }
    else if(edge_weight_type=="EUC_2D" || edge_weight_type=="CEIL_2D" || edge_weight_type=="ATT" || edge_weight_type=="GEO"){
        
        _DM=(DistanceMode)DistanceMode::Euclidian;
        if(edge_weight_type=="CEIL_2D"){
            _DM=(DistanceMode)DistanceMode::Pseudoeuclidian;
        }
        else if(edge_weight_type=="ATT"){
            _DM=(DistanceMode)DistanceMode::Fischetti;
        }
        else if(edge_weight_type=="GEO"){
            _DM=(DistanceMode)DistanceMode::Geographical;
        }
    }
    else{
        throw runtime_error("Edge weight type "+edge_weight_type+" is not implemented");
    }
    
}

void Instance::readExplicit(ifstream& input, string edge_weight_format)
{
    if(!gotoSection(input, "EDGE_WEIGHT_SECTION")){
        throw runtime_error("Could not find EDGE_WEIGHT_SECTION in file "+ instancename());
    }

    _D=new double*[_nnodes];
    for(int i=0;i<_nnodes;i++){
        _D[i]=new double[_nnodes];
        for(int j=0;j<i;j++){
            _D[i][j]=0;
        }
    }

    if(edge_weight_format=="FULL_MATRIX"){
        for(int i=0;i<_nnodes;i++){
            for(int j=0;j<_nnodes;j++){
                double d;
                input>>d;
                setTravellingTime(i,j,d);
            }
        }
    }
    else if (edge_weight_format=="UPPER_ROW"){
        for(int i=0;i<_nnodes-1;i++){
            for(int j=i+1;j<_nnodes;j++){
                double d;
                input>>d;
                setTravellingTime(i,j,d);
            }
        }
    }
    else if(edge_weight_format=="LOWER_ROW"){
        for(int i=1;i<_nnodes;i++){
            for(int j=0;j<i;j++){
                double d;
                input>>d;
                setTravellingTime(i,j,d);
            }
        }
    }
    else if(edge_weight_format=="UPPER_DIAG_ROW"){
        for(int i=0;i<_nnodes;i++){
            for(int j=i;j<_nnodes;j++){
                double d;
                input>>d;
                setTravellingTime(i,j,d);
            }
        }
    }
    else if(edge_weight_format=="LOWER_DIAG_ROW"){
        for(int i=0;i<_nnodes;i++){
            for(int j=0;j<=i;j++){
                string d;
                input>>d;
                cout << d << " ";
                setTravellingTime(i,j,stod(d));
            }
        }
        
    }
    else{
        throw runtime_error("EDGE_WEIGHT_FORMAT "+edge_weight_format+" is not implemented");
    }
}


Instance::~Instance()
{
    for(int i=0;i<_nnodes;i++) delete _Nodes[i];
    delete[] _Nodes;
    if(_D!=NULL){
        for(int i=1;i<_nnodes;i++) delete[] _D[i];
        delete[] _D;
    }
}

void Instance::print()
{
	//TODO
    cout<<"File: "<<_name<<endl;
    cout<<"DM= "<<(int)_DM<<endl;
    for(int i=0;i<_nnodes;i++) _Nodes[i]->print();
}

double Instance::travellingtime(node* Ni, node* Nj){
    double dij = 0;
    if (_DM == DistanceMode::Euclidian)
    {
        dij= round(sqrt(pow(Ni->_x - Nj->_x, 2) + pow(Ni->_y - Nj->_y, 2)));
    }
    else if (_DM == DistanceMode::Pseudoeuclidian)
    {
        
        double rij = sqrt((pow(Ni->_x - Nj->_x, 2) + pow(Ni->_y - Nj->_y, 2))/10);
        //int tij = int(rij);//floor
        //if (tij < rij) dij = tij + 1;
        //else dij = tij;
        dij = ceil(rij); //It was a "floor" -> changed to "ceil" (seems consistent with OP)
    }
    else if (_DM == DistanceMode::Fischetti)
    {
        double rij = 100*sqrt((pow(Ni->_x - Nj->_x, 2) + pow(Ni->_y - Nj->_y, 2)));
        dij= (rij);
        
    }
    else if (_DM == DistanceMode::Geographical)
    {
        //Convert to lat and lng
        
        double lat_i = coordtogeo(Ni->_x);
        double lng_i = coordtogeo(Ni->_y);

        double lat_j = coordtogeo(Nj->_x);
        double lng_j = coordtogeo(Nj->_y);

        double rrr = 6378.388; //Earth radius
        double q1 = cos(lng_i  - lng_j );
        double q2 = cos(lat_i  - lat_j );
        double q3 = cos(lat_i  + lat_j );
        dij = (int)(rrr * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }
    else if(_DM == DistanceMode::Explicit)
    {
        if(Ni->id()==Nj->id()) return 0;
        int i = max(Ni->id(), Nj->id());
        int j = min(Ni->id(), Nj->id());
        dij=_D[i-1][j-1];
    }
    else
    {
        throw runtime_error("Error: Distance Mode is not implemented");
    }
    return dij;
}

//double Instance::travellingtimeedge() 
//{
//    double de = 0;
//    return de;
//}