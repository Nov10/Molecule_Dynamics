using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Linq;


class Bond
{
    public int Particle1;
    public int Particle2;
    public float Distance;
    public float Strength;

    public Bond(int particle1, int particle2, float distance, float stength)
    {
        Particle1 = particle1;
        Particle2 = particle2;
        Distance = distance;
        Strength = stength;
    }
}

class Angle
{
    public int Particle1;
    public int Particle2;
    public int Particle3;
    public float Angle0;
    public float KAngle;

    public Angle(int particle1, int particle2, int particle3, float angle0, float kAngle)
    {
        Particle1 = particle1;
        Particle2 = particle2;
        Particle3 = particle3;
        Angle0 = angle0;
        KAngle = kAngle;
    }
}



public class SimulationSolver : MonoBehaviour
{ 
    Molecule[] Molecules;
    public int Numbers;
    public float T;
    List<Bond> Bonds;
    List<Angle> Angles;
    public float[] mols;
    public float[] chrg;
    [SerializeField] Transform Parent;
    [SerializeField] Transform Sphere1;
    [SerializeField] Transform Sphere2;
    Transform[] Visualizers;
    //Constants
    const float kb = 0.8314459920816467f;
    const float NA = 6.0221409e+26f;
    const float ech = 1.60217662E-19f;
    const float kc = 8.9875517923E9f * NA * 1E30f * ech * ech / 1E24f;


    [SerializeField] Vector3 BoundaryCenter;
    [SerializeField] float BoundaryHalfSize;

    

    private void Awake()
    {
        Molecules = new Molecule[Numbers];
        Visualizers = new Transform[Numbers];
        for(int i = 0; i < Numbers; i++)
        {
            Molecule molecule = new Molecule();
            molecule.InitializeRandom();

            Molecules[i] = molecule;
            //Visualizers[i] = Instantiate(Sphere, Parent);
        }

        Bonds = new List<Bond>();
        for(int j = 0; j < Numbers/3; j++)
        {
            Visualizers[3*j] = Instantiate(Sphere1, Parent);
            Visualizers[3*j+2] = Instantiate(Sphere1, Parent);
            Visualizers[3*j+1] = Instantiate(Sphere2, Parent);
            Bonds.Add(new Bond(3 * j, 3 * j + 1, 1.0f, 100.0f));
            Bonds.Add(new Bond(3 * j + 1, 3 * j + 2, 1.0f, 100.0f));
        }

        Angles = new List<Angle>();
        for(int i = 0; i < Numbers/3; i++)
        {
            Angles.Add(new Angle(3 * i, 3 * i + 1, 3 * i + 2, 109.47f * Mathf.PI/180.0f, 35300.0f / 2));
        }

        mols = new float[Numbers];
        for (int i = 0; i < Numbers / 3; i++)
        {
            mols[3 * i] = i;
            mols[3 * i + 1] = i;
            mols[3 * i + 2] = i;
        }

        chrg = new float[Numbers];
        for (int i = 0; i < Numbers / 3; i++)
        {
            chrg[3 * i] = 0.41f;
            chrg[3 * i + 1] = -2.82f;
            chrg[3 * i + 2] = 0.41f;
        }


    }

    private void Update()
    {
        float deltaTime = Time.deltaTime;

        float sum_Ek = 0;
        float avEK = 0;

        LoopForMoleculels((Molecule molecule) =>
        {
            molecule.Update(deltaTime, Molecules, BoundaryCenter, Vector3.one * BoundaryHalfSize);
            sum_Ek += molecule.CalculateEk();
        });

        avEK = sum_Ek / Numbers;

        LoopForMoleculels((Molecule molecule) =>
        {
            molecule.totalEk = avEK;
            molecule.Velocity = molecule.Velocity.normalized * (float)molecule.rescaleT(molecule.Velocity.magnitude, T);
        });

        Vector3[] dBE = dBEpot();
        Vector3[] dba = dBA();
        LoopForMoleculels((int i) =>
        {
            Molecules[i].Velocity -= dBE[i].sqrMagnitude <= 1 ? dBE[i] : dBE[i].normalized;
            Molecules[i].Velocity -= dba[i].sqrMagnitude <= 1 ? dba[i] : dba[i].normalized;
            //Molecules[i].Velocity -= dba[i];
        });

        LoopForMoleculels((int i) =>
        {
            Vector3 v = (coul(i) / Molecules[i].m) * deltaTime;
            Molecules[i].Velocity += v.sqrMagnitude <= 1 ? v : v.normalized;

        });


        Visualize();
    }

    private void OnDrawGizmos()
    {
        Gizmos.color = Color.white;
        Gizmos.DrawWireCube(BoundaryCenter, Vector3.one * BoundaryHalfSize * 2);
    }

    void Visualize()
    {
        LoopForMoleculels((int idx) =>
        {
            Visualizers[idx].position = Molecules[idx].Position;
        });

        for (int i = 0; i < Numbers; i++)
        {
            for (int j = 0; j < Bonds.Count; j++)
            {
                if (Bonds[j].Particle1 == i || Bonds[j].Particle2 == i)
                {
                    int ii = 0;
                    if (Bonds[j].Particle1 == i)
                    {
                        ii = Bonds[j].Particle2;
                    }
                    else
                    {
                        ii = Bonds[j].Particle1;
                    }
                    Debug.DrawLine(Molecules[i].Position, Molecules[ii].Position, Color.red);
                }
            }
        }
    }
    void LoopForMoleculels(System.Action<int> function)
    {
        for(int i = 0; i < Numbers; i++)
        {
            function(i);
        }
    }

    void LoopForMoleculels(System.Action<Molecule> function)
    {
        for (int i = 0; i < Numbers; i++)
        {
            function(Molecules[i]);
        }
    }


    float[] Bepot()
    {
        float[] bps = new float[Numbers];
        for(int i = 0; i < Numbers; i++)
        {
            for(int j = 0; j < Bonds.Count; j++)
            {
                if(Bonds[j].Particle1 == i || Bonds[j].Particle2 == i)
                {
                    int ii = 0;
                    if(Bonds[j].Particle1 == i)
                    {
                        ii = Bonds[j].Particle2;
                    }
                    else
                    {
                        ii = Bonds[j].Particle1;
                    }
                    float dr0 = Bonds[j].Distance;
                    float e0 = Bonds[j].Strength;
                    float dr = (Molecules[i].Position - Molecules[ii].Position).magnitude;
                    float BE = e0 * (dr - dr0) * (dr - dr0);
                    bps[i] += 0.5f * BE;
                }
            }
        }
        return bps;
    }

    Vector3[] dBEpot()
    {
        Vector3[] bps = new Vector3[Numbers];
        for (int i = 0; i < Numbers; i++)
        {
            for (int j = 0; j < Bonds.Count; j++)
            {
                if (Bonds[j].Particle1 == i || Bonds[j].Particle2 == i)
                {
                    int ii = 0;
                    if (Bonds[j].Particle1 == i)
                    {
                        ii = Bonds[j].Particle2;
                    }
                    else
                    {
                        ii = Bonds[j].Particle1;
                    }
                    float dr0 = Bonds[j].Distance;
                    float e0 = Bonds[j].Strength;
                    float dr = (Molecules[i].Position - Molecules[ii].Position).magnitude;
                    Vector3 dBE = 2* e0 * (dr - dr0) * (Molecules[i].Position - Molecules[ii].Position)/dr;
                    bps[i] += 0.5f * dBE;
                }
            }
        }
        return bps;
    }

    Vector3[] dBA()
    {
        Vector3[] aps = new Vector3[Numbers];
        int a1 = 0;
        int a2 = 0;
        int a3 = 0;
        float th00 = 0.0f;
        float e0 = 0.0f;
        Vector3 r1;
        Vector3 r2;
        float ar1;
        float ar2;
        float dot;
        float ndot;
        float th;
        float dUdth;
        Vector3 numerator;
        float denominator;
        Vector3 dUdr;
        Vector3 n1;
        Vector3 n2;
        Vector3 n3;

        for(int i = 0; i < Numbers; i++)
        {
            for (int j = 0; j < Angles.Count; j++)
            {
                Debug.Log("A");
                a1 = Angles[j].Particle1;
                a2 = Angles[j].Particle2;
                a3 = Angles[j].Particle3;
                if (i == a1 || i == a2 || i == a3)
                {
                    th00 = Angles[j].Angle0;
                    e0 = Angles[j].KAngle;

                }
                if (i == a1 || i == a2)
                {
                    r1 = Molecules[a1].Position - Molecules[a2].Position;
                    r2 = Molecules[a3].Position - Molecules[a2].Position;
                    ar1 = (Molecules[a1].Position - Molecules[a2].Position).magnitude;
                    ar2 = (Molecules[a3].Position - Molecules[a2].Position).magnitude;
                }
                else
                {
                    r1 = Molecules[a3].Position - Molecules[a2].Position;
                    r2 = Molecules[a1].Position - Molecules[a2].Position;
                    ar1 = (Molecules[a3].Position - Molecules[a2].Position).magnitude;
                    ar2 = (Molecules[a1].Position - Molecules[a2].Position).magnitude;
                }
                dot = r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2];
                ndot = dot / (ar1 * ar2);
                th = Mathf.Acos(ndot);
                dUdth = -2.0f * e0 * (th - th00);
                if (a1 == i || a3 == i)
                {
                    numerator = (r2 / (ar1 * ar2)) - Vector3.one * (dot / (ar1 * ar1 * ar1 * ar2 * 2.0f));
                    denominator = Mathf.Sqrt((1.0f - ndot * ndot));
                    dUdr = dUdth * numerator / denominator;
                    aps[i] += dUdr;
                }
                if (i == a2)
                {
                    denominator = Mathf.Sqrt(1.0f - ndot * ndot);
                    n1 = -(r2 + r1);
                    n2 = dot * r1 / (ar1 * ar1);
                    n3 = dot * r2 / (ar2 * ar2);
                    numerator = (n1 + n2 + n3) / (ar1 * ar2);
                    dUdr = dUdth * numerator / denominator;
                    aps[i] += dUdr;
                }
            }
        }
        return aps;
    }

    Vector3 coul(int i)
    {
        float q0 = chrg[i];
        float[] qs = chrg;
        for (int j = 0; j < Numbers; j++)
        {
            if (mols[i] == mols[j])
            {
                qs[j] = 0.0f;
            }
        }
        List<float> ListQs = qs.ToList();
        ListQs.RemoveAt(i);
        qs = ListQs.ToArray();
        List<Vector3> deltaPositions = new List<Vector3>();
        LoopForMoleculels((int m) =>
        {
        if (i != m)
            {
                deltaPositions.Add(Molecules[m].Position - Molecules[i].Position);
            }
        });

        Vector3 F = Vector3.zero;

        for(int k = 0; k < deltaPositions.Count; k++)
        {
            F += ((q0 * ListQs[k] * kc) / deltaPositions[k].sqrMagnitude) * deltaPositions[k];
        }
        return F;
    }



}
