
// This is a C++ file with a simple model used to test the correct installation of CPLEX.
// It defines a basic mixed-integer linear programming (MILP) problem and solves it using CPLEX.
#include <ilcplex/ilocplex.h>

ILOSTLBEGIN


int main (void) {
   IloEnv env;
   try {
      IloModel model(env);

      IloNumVarArray var(env);
      IloRangeArray con(env);
      
      var.add(IloNumVar(env, 0.0, 40.0));
      var.add(IloNumVar(env));
      var.add(IloNumVar(env));
      var.add(IloNumVar(env, 2.0, 3.0, ILOINT));
      model.add(IloMaximize(env, var[0] + 2 * var[1] + 3 * var[2] + var[3]));

      con.add(-var[0] + var[1] + var[2] + 10 * var[3] <= 20);
      con.add(var[0] - 3 * var[1] + var[2] <= 30);
      con.add(var[1] - 3.5 * var[3] == 0);
      model.add(con);

      IloCplex cplex(model);
      cplex.solve();

      env.out() << "Solution status = " << cplex.getStatus() << endl;
      env.out() << "Solution value  = " << cplex.getObjValue() << endl;

      IloNumArray vals(env);
      cplex.getValues(vals, var);
      env.out() << "Values        = " << vals << endl;
      cplex.getSlacks(vals, con);
      env.out() << "Slacks        = " << vals << endl;

	  cplex.exportModel("milpsample.lp");
   }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }

   env.end();
   return 0;

}


