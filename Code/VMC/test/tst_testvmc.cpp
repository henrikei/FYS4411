#include <QtCore/QString>
#include <QtTest/QtTest>
#include "waveFunction/helium.h"
#include <armadillo>

using namespace arma;

class TestVMC : public QObject
{
    Q_OBJECT
    
public:
    TestVMC();
    
private Q_SLOTS:
    void initTestCase();
    void cleanupTestCase();
    void testwavefunction();
};

TestVMC::TestVMC()
{
}

void TestVMC::initTestCase()
{
}

void TestVMC::cleanupTestCase()
{
}

void TestVMC::testwavefunction()
{
    helium a;
    a.setAlpha(0.9);
    mat pos = zeros(2,3);
    cout << a.getValue(pos) << endl;
    QVERIFY2(a.getValue(pos)==1, "Failure");

    pos << 1 << 2 << 3 << endr
        << 2 << 3 << 1;
    cout << pos << endl;
    cout << a.getValue(pos) << endl;
    QVERIFY2(abs(a.getValue(pos) - 0.0011886) < 1e-5 , "Failure");
}

QTEST_APPLESS_MAIN(TestVMC)

#include "tst_testvmc.moc"
