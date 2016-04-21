#include "MachineLearning.h"
#include "MachineLearningOthers.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;

template<typename DATA> void readWineData(DATA& result)
{
    ifstream fin("Datasets/wine.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        --label;
        if(label < 0 || label > 2) continue;
        data >> comma;
        for(int i = 0; i < 12; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readStatlogData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/sat.trn" : "Datasets/sat.tst");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 36; i++)
        {
            double currentDigit;
            data >> currentDigit;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label == 7) label = 6;
        if(label < 0 || label > 6) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readSpectData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/SPECT.train" : "Datasets/SPECT.test");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        data >> comma;
        for(int i = 0; i < 22; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readSpamData(DATA& result)
{
    ifstream fin("Datasets/spambase.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 57; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label != 0 && label != 1) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readPimaData(DATA& result)
{
    ifstream fin("Datasets/pima-indians-diabetes.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 8; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label != 0 && label != 1) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readLetterData(DATA& result)
{
    ifstream fin("Datasets/letter-recognition.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        char sLabel;
        data >> sLabel;
        int label = sLabel - 'A';
        if(label < 0 || label > 26) continue;
        data >> comma;
        for(int i = 0; i < 16; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readIonosphereData(DATA& result)
{
    ifstream fin("Datasets/ionosphere.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 34; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        string sLabel;
        data >> sLabel;
        int label;
        if(sLabel == "g") label = 0;
        else if(sLabel == "b") label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readGlassData(DATA& result)
{
    ifstream fin("Datasets/glass.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 10; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label > 7 || label < 1) continue;
        result.addZ(x, label - 1);
    }
}

template<typename DATA> void readMadelonData(DATA& result, bool isTrain)
{;
    ifstream fin(isTrain ? "Datasets/madelon_train.data" : "Datasets/madelon_valid.data");
    ifstream fin2(isTrain ? "Datasets/madelon_train.labels" : "Datasets/madelon_valid.labels");
    while(!fin.eof() && !fin2.eof())
    {
        string line, line2;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 500; i++)
        {
            double currentDigit;
            data >> currentDigit;

            x.append( currentDigit );
        }
        assert(x.getSize() == 500);
        getline(fin2, line2, '\n');
        stringstream data2;
        data2 << line2;
        int dummy, label;
        data2 >> dummy;
        if(dummy == -1) label = 0;
        else if(dummy == 1) label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readArceneData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/arcene_train.data" : "Datasets/arcene_valid.data");
    ifstream fin2(isTrain ? "Datasets/arcene_train.labels" : "Datasets/arcene_valid.labels");
    while(!fin.eof() && !fin2.eof())
    {
        string line, line2;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 10000; i++)
        {
            double currentDigit;
            data >> currentDigit;

            x.append( currentDigit );
        }
        assert(x.getSize() == 10000);
        getline(fin2, line2, '\n');
        stringstream data2;
        data2 << line2;
        int dummy, label;
        data2 >> dummy;
        if(dummy == -1) label = 0;
        else if(dummy == 1) label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readWDBCData(DATA& result)
{
    ifstream fin("Datasets/wdbc.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int dummy;
        data >> dummy;
        data >> comma;
        char sLabel;
        data >> sLabel;
        int label;
        if(sLabel == 'B') label = 0;
        else if(sLabel == 'M') label = 1;
        else continue;
        data >> comma;
        for(int i = 0; i < 30; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        assert(x.getSize() == 30);
        result.addZ(x, label);
    }
}

template<typename DATA> void readBanknoteData(DATA& result)
{
    ifstream fin("Datasets/data_banknote_authentication.txt");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 4; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        assert(x.getSize() == 4);
        result.addZ(x, label);
    }
}

template<typename DATA> void readCNEAData(DATA& result)
{
    ifstream fin("Datasets/CNAE-9.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        data >> comma;
        for(int i = 0; i < 856; i++)
        {
            double currentDigit;
            data >> currentDigit;
            data >> comma;
            x.append( currentDigit);
        }
        assert(x.getSize() == 856);
        result.addZ(x, label - 1);
    }
}

template<typename DATA> void readIrisData(DATA& result)
{
    ifstream fin("Datasets/iris.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 4; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        string sLabel;
        data >> sLabel;
        int label;
        if(sLabel == "Iris-setosa") label = 0;
        else if(sLabel == "Iris-versicolor") label = 1;
        else if(sLabel == "Iris-virginica") label = 2;
        else continue;
        result.addZ(x, label);
    }
}


template<typename DATA> void readDigitData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain? "optdigits.tra" : "optdigits.tes");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 64; i++)
        {
            double currentDigit;
            data >> currentDigit;
            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        result.addZ(x, label);
    }
}

Matrix<double> sampleCostDeterministic(int k)
{
    assert(k > 0);
    int count = 0;
    Matrix<double> result(k, k);
    for(int r = 0; r < k; ++r)
        for(int c = 0; c < k; ++c) if(r != c) result(r, c) = count++ % 2 ? 0.01 : 1;
    scaleCostMatrix(result);
    return result;
}

template<typename T> pair<PermutedData<T>, PermutedData<T> >
makeData(Vector<T>& dataM)
{
    return createTrainingTestSetsStatified(dataM.lastItem());
}

template<typename T> pair<PermutedData<T>, PermutedData<T> >
makeDataDivided(Vector<T>& dataM)
{
    PermutedData<T> p1(dataM[dataM.getSize() - 2]);
    for(int i = 0; i < dataM[dataM.getSize() - 2].getSize(); ++i) p1.addIndex(i);
    PermutedData<T> p2(dataM.lastItem());
    for(int i = 0; i < dataM.lastItem().getSize(); ++i) p2.addIndex(i);
    return make_pair(p1, p2);
}

template<typename LEARNER> int testNumericalClassifier()
{
    DEBUG("Started Reading");
    typedef InMemoryData<NUMERIC_X, int> T;
    Vector<T> dataM(50);//make many enough to avoid ref realloc
    Vector<pair<PermutedData<T>, PermutedData<T> > > data;

    dataM.append(T());
    readIrisData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    /*dataM.append(T());
    readDigitData(dataM.lastItem(), true);
    dataM.append(T());
    readDigitData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));*/

    /*dataM.append(T());
    readCNEAData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readBanknoteData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readWDBCData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readGlassData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readIonosphereData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    dataM.append(T());
    readLetterData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readPimaData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readSpamData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readSpectData(dataM.lastItem(), true);
    dataM.append(T());
    readSpectData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readStatlogData(dataM.lastItem(), true);
    dataM.append(T());
    readStatlogData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readWineData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    dataM.append(T());
    readArceneData(dataM.lastItem(), true);
    dataM.append(T());
    readArceneData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readMadelonData(dataM.lastItem(), true);
    dataM.append(T());
    readMadelonData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));*/

    DEBUG("Done Reading");
    int reportNumber = time(0);
    string fAcc = "reportAcc" + toString(reportNumber) + ".csv";
    string fAveAcc = "reportAveAcc" + toString(reportNumber) + ".csv";
    string fTimer = "reportTimer" + toString(reportNumber) + ".csv";
    string fCount = "reportFeature" + toString(reportNumber) + ".csv";
    string fCost = "reportCost" + toString(reportNumber) + ".csv";
    ++reportNumber;
    for(int i = 0; i < data.getSize(); ++i)
    {
        if(false)//cost
        {
            /*int start = clock();
            int k = findNClasses(data[i].first);
            Matrix<double> c = sampleCostDeterministic(k);
            for(int z = 0; z < k; ++z)
            {
                for(int j = 0; j < k; ++j)
                    cout << c(z, j) << " ";
                cout << endl;
            }
            //LEARNER s(T(data[i].first).data);
            LEARNER s(T(data[i].first).data, c);
            Matrix<int> confusion = evaluateConfusion(evaluateLearner<int>(s, T(data[i].second).data));
            for(int z = 0; z < k; ++z)
            {
                for(int j = 0; j < k; ++j)
                    cout << confusion(z, j) << " ";
                cout << endl;
            }
            ClassifierStats cs(confusion);
            double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
            double cost = evalConfusionCost(confusion, c);
            DEBUG(cost);
            cs.debug();
            addToCSV(Vector<string>(1, toString(cost)), fCost.c_str());
            addToCSV(Vector<string>(1, toString(timediff)), fTimer.c_str());*/
        }
        else
        {
            int start = clock();
            LEARNER NeuralNetworkIris(data[i].first);
            ClassifierStats cs(evaluateConfusion(evaluateLearner<int>(NeuralNetworkIris, data[i].second)));
            double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
            cs.debug();

            /*addToCSV(Vector<string>(1, toString(cs.acc.mean)), fAcc.c_str());
            addToCSV(Vector<string>(1, toString(cs.bac.mean)), fAveAcc.c_str());
            addToCSV(Vector<string>(1, toString(timediff)), fTimer.c_str());*/
            //addToCSV(Vector<string>(1, toString(s.model.f.fMap.getSize())), fCount.c_str());
        }
        //system("PAUSE");
    }
    return 0;
}

void testNumericalClassifiers()
{
//    testNumericalClassifier<SSVM>();
//    testNumericalClassifier<DecisionTree>();
//
//
//    testNumericalClassifier<SNN>();
//    testNumericalClassifier<SLSVM>();
//    testNumericalClassifier<RMBoost<> >();
//    testNumericalClassifier<AdaBoostSamme<> >();
//    testNumericalClassifier<RandomForest>();
//
//
//    testNumericalClassifier<SImbSVM>();
//
//    testNumericalClassifier<SRaceLSVM>();
//    testNumericalClassifier<SOnlineNN>();
//
//    testNumericalClassifier<ScaledLearner<NoParamsLearner<KNNClassifier<>, int>, int> >();
//
//
//
    testNumericalClassifier<SimpleBestCombiner>();

    //inactive learners
    //testNumericalClassifier<ScaledLearner<MeanNN<>, int> >();
    //testNumericalClassifier<NumericalBayes>();

    //feature selection
    //testNumericalClassifier<SmartFSLearner<> >();

    //cost learning
    //testNumericalClassifier<SBoostedCostSVM>();
    //testNumericalClassifier<SAveCostSVM>();
    //testNumericalClassifier<CostLearner<> >();

}





























//holds the digit data
struct Digit
{
    vector< vector < int > > color;
    int actualDigit;
    Digit(){}
    Digit(vector< vector < int > > color_, int actualDigit_)
    :color(color_),actualDigit(actualDigit_)
    {}
};

//reads in the digits
vector< Digit > readDigits( string filename)
{
    //the data format is 65 comma-separated integer values per line
    //the last value is the actual digit, the 64 are the colors

    vector< Digit > result;
    ifstream fin(filename.c_str());
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        Digit readDigit;
        for(int i = 0; i < 8; i++)
        {
            vector< int > row;
            for(int j = 0; j < 8; j++)
            {
                int currentDigit;
                data >> currentDigit;
                char comma;
                data >> comma;
                row.push_back( currentDigit );
            }
            readDigit.color.push_back(row);
        }
        data >> readDigit.actualDigit;
        result.push_back(readDigit);
    }


    return result;
}

void testKMeans()
{
    Vector<Point2> points;
    int N = 20000;
    for(int i = 0; i < N; ++i)
    {
        points.append(Point2(GlobalRNG.uniform01(), GlobalRNG.uniform01()));
    }
    Vector<int> result = KMeans<Point2>::findClusters(points, 10);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i]);
    }
}

void DDDKMeans()
{
    InMemoryData<NUMERIC_X, int> data;
    readIrisData(data);
    Vector<NUMERIC_X> points;
    for(int i = 0; i < data.getSize(); ++i)
    {
        points.append(data.getX(i));
    }
    Vector<int> result = KMeans<NUMERIC_X>::findClusters(points, 3);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i]);
    }
}

void testKMeans2()
{
    vector< Digit > testSet = readDigits( "optdigits.tra" );
    Vector<Point<double, 64> > points;
    for(int k = 0; k < testSet.size(); k++)//
    {
        Point<double, 64> instance;
        for(int i = 0; i < 8; i++)
            for(int j = 0; j < 8; j++)
                instance[i * 8 + j] = testSet[k].color[i][j];
        points.append(instance);
    }
    Vector<int> result = KMeans<Point<double, 64> >::findClusters(points, 10);
    for(int i = 0; i < result.getSize(); ++i)
    {
        DEBUG(result[i]);
    }
}

void testAPriori()
{
    Vector<Vector<int> > baskets;
    Vector<int> b1, b2, b3, b4;
    b1.append(0);
    b1.append(1);
    b1.append(2);
    b1.append(3);
    baskets.append(b1);
    b2.append(5);
    b2.append(1);
    b2.append(2);
    b2.append(4);
    baskets.append(b2);
    b3.append(7);
    b3.append(1);
    b3.append(2);
    b3.append(6);
    baskets.append(b3);
    b4.append(1);
    b4.append(0);
    b4.append(4);
    b4.append(6);
    baskets.append(b4);
    APriori ap;
    ap.noCutProcess(baskets, 3);
    for(LcpTreap<Vector<int>, int>::Iterator i(ap.counts.begin()); i != ap.counts.end(); ++i)
    {
        for(int j = 0; j < i->key.getSize(); ++j)
        {
            DEBUG(i->key[j]);
        }
        DEBUG(i->value);
    }
}


struct GridWorld
{
    DiscreteValueFunction u;
    int state, nEpisodes, nextState;
    double reward()
    {
        if(state == 3) return 1;
        if(state == 7) return -1;
        return -0.04;
    }
    double discountRate(){return 1;}
    double goToNextState(){state = nextState;}
    double pickNextState()
    {
        int row = state % 4, column = state / 4;
        int rows[4] = {row+1,row,row,row-1};
        int columns[4] = {column,column+1,column-1,column};
        bool set = false;
        for(int i = 0; i < 4; ++i)
        {
            if(rows[i] >= 0 && rows[i] <= 3 && columns[i] >= 0 && columns[i] <= 2 && !(rows[i] == 1 && columns[i] == 1))
            {
                int newState = rows[i] + columns[i] * 4;
                assert(state !=newState);
                if(!set || u.values[newState].first > u.values[nextState].first) {nextState = newState; set = true;}
            }
        }
        assert(set);
        return u.values[nextState].first;
    }
    bool isInFinalState(){return state == 3 || state == 7;}
    double learningRate(){return u.learningRate(state);}
    bool hasMoreEpisodes(){return nEpisodes;}
    double startEpisode()
    {
        do{state = GlobalRNG.mod(12);} while(state == 5);
        --nEpisodes;
        return u.values[state].first;}
    void updateCurrentStateValue(double delta){u.updateValue(state, delta);}
    GridWorld():nEpisodes(100), u(12){}
    void debug()
    {
        for(int i = 0; i < 3; ++i)
        {
            for(int j = 0; j < 4; ++j)
            {
                cout << " " << u.values[j + i * 4].first;
            }
            cout << endl;
        }
    }
};

void testReinforcement()
{
    GridWorld g;
    TDLearning(g);
    g.debug();
}

int main(int argc, char *argv[])
{
    DDDKMeans();
    //testKMeans2();
    //testReinforcement();
    //testAPriori();
    //for(int i = 0; i < 1; ++i) testNumericalClassifiers();
	return 0;
}


