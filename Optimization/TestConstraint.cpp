#include "ContraintProcessing.h"
#include "../Utils/Debug.h"
using namespace igmdk;
using namespace std;

void testAC3()
{
    //Example: trace execution on (0)<->(0,1)<->(0,1,2)<->(0,1,2,3) with <->
    //denoting alldiff constraint
    AllDifferent ad;

    ConstraintGraph<AllDifferent::Handle> cg;
    for(int i = 0; i < 4; ++i)
    {
        cg.addVariable(i+1);
        cg.variables[i].setAll(true);
        ad.addVariable(i);
    }
    for(int i = 0; i < 4; ++i) cg.variables[i].output();
    for(int i = 0; i < 3; ++i) cg.addConstraint(i, i+1, ad.handle);
        //for(int j = i + 1; j < 4; ++j)
        //    cg.addConstraint(i, j, AllDifferent());
    DEBUG(cg.AC3());
    for(int i = 0; i < 4; ++i) cg.variables[i].output();
    //Vector<int> result;
    //backtrackFind(cg.variables, result, AllDifferent());
    //for(int i = 0; i < result.getSize(); ++i) DEBUG(result[i]);
}

struct Sudoku
{
    AllDifferent ad[3][9];
    ConstraintGraph<AllDifferent::Handle> cg;
    Sudoku(int* values)
    {
        for(int i = 0; i < 81; ++i)
        {
            cg.addVariable(9);
            if(values[i])
            {
                cg.variables[i].setAll(false);
                cg.variables[i].set(values[i]-1, true);
            }
            else cg.variables[i].setAll(true);
        }
        for(int i = 0; i < 9; ++i)
        {
            int rowStart = i * 9, columnStart = i, boxStart = (i / 3) * 27 + (i % 3) * 3;
            for(int j = 0; j < 9; ++j)
            {
                int rowMember = rowStart + j;
                int columnMember = columnStart + j*9;
                int boxMember = boxStart + (j / 3) * 9 + j % 3;
                ad[0][i].addVariable(rowMember);
                ad[1][i].addVariable(columnMember);
                ad[2][i].addVariable(boxMember);
                if(j == 8) continue;
                for(int k = j+1; k < 9; ++k)
                {
                    int boxMember2 = boxStart + (k / 3) * 9 + k % 3;
                    cg.addConstraint(rowMember, rowStart + k, ad[0][i].handle);
                    cg.addConstraint(columnMember, columnStart + k*9, ad[1][i].handle);
                    cg.addConstraint(boxMember, boxMember2, ad[2][i].handle);
                }
            }
        }
        DEBUG(cg.AC3());
        for(int i = 0; i < 81; ++i)
        {
            if(i % 9 == 0) cout << endl;
            int count = 0, value = -1;
            for(int j = 0; j < 9; ++j)
            {
                if(cg.variables[i][j])
                {
                    ++count;
                    value = j + 1;
                }
            }
            if(count > 1) cout << "x";
            else cout << value;

        }
        cout << endl;
    }
};

void testSudoku()
{
    int hard[] =
    {
        0,1,2,0,6,0,8,0,0,
        0,0,0,0,3,0,0,5,0,
        6,0,0,4,0,0,0,0,7,
        0,0,6,0,0,0,0,1,0,
        0,9,7,0,0,0,6,4,0,
        0,8,0,0,0,0,7,0,0,
        8,0,0,0,0,1,0,0,3,
        0,4,0,0,5,0,0,0,0,
        0,0,1,0,2,0,9,7,0
    };
    int array[] =
    {
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0
    };
    int easy[] =
    {
        2,0,1,7,9,5,0,0,0,
        0,0,9,0,0,8,0,1,0,
        0,0,0,3,0,1,0,0,7,
        0,2,0,5,0,0,1,7,8,
        0,8,0,0,0,0,0,9,0,
        1,5,7,0,0,4,0,3,0,
        6,0,0,8,0,2,0,0,0,
        0,9,0,6,0,0,5,0,0,
        0,0,0,1,7,9,3,0,6
    };
    int medium[] =
    {
        0,0,5,0,0,3,2,9,0,
        9,0,0,2,0,0,0,3,4,
        0,0,0,0,1,0,8,0,0,
        0,0,0,0,9,0,0,7,1,
        0,0,0,6,0,5,0,0,0,
        7,3,0,0,2,0,0,0,0,
        0,0,7,0,6,0,0,0,0,
        6,8,0,0,0,9,0,0,2,
        0,5,2,8,0,0,6,0,0
    };
    Sudoku sd(easy);
}


int main()
{
    testSudoku();
    testAC3();
	return 0;
}
