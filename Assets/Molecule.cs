using System;
using Unity.VisualScripting;
using UnityEngine;

[System.Serializable]
public class Molecule
{
    public Vector3 Position;
    public Vector3 Velocity;
    public float Radius = 1;
    public float m; //mass of paarticles in Daltons
    public float T0; //temperature in K
    public float k0;  //Boltzmann constant
    public float k1; //molar gas constant
    public float k2;
    public float k3;
    public float kb;
    public float eps = 1;
    public float sig = 1;
    public float totalEk;

    public Molecule()
    {
        Position = new Vector3(UnityEngine.Random.Range(0, 1f), UnityEngine.Random.Range(0, 1f), UnityEngine.Random.Range(0, 1f));
        Velocity = new Vector3(UnityEngine.Random.Range(0, 1f), UnityEngine.Random.Range(0, 1f), UnityEngine.Random.Range(0, 1f));
        k0 = 1.38064852E-23f;
        T0 = 10;
        m = 1;
        k1 = k0 * 6.0221409e+26f;
        k2 = k1 * 1.0E20f;
        k3 = k2 / 1.0E24f;
        kb = k3;
    }


    public void InitializeRandom()
    {
        Position = new Vector3(UnityEngine.Random.Range(-1, 1f), UnityEngine.Random.Range(-1, 1f), UnityEngine.Random.Range(-1, 1f)) * 10f;
        Velocity = new Vector3(0, 0, 0);
    }

    public float CalculateEk()
    {
        return 0.5f * m * Velocity.sqrMagnitude;
    }

    public float rescaleT(float v, float T)
    {
        float Tnow;
        Tnow = (2.0f / 3) * totalEk / kb;
        float lam = Mathf.Sqrt(T / Tnow);
        float vnew = lam * v;
        return vnew;
    }

    float CalculatePotentialEnergy(Molecule[] moles)
    {
        float sum = 0;
        for(int i = 0; i < moles.Length; i++)
        {
            if (moles[i] != this)
            {
                sum += CalculatePursion(moles[i]);
            }
        }
        return sum;
    }

    Vector3 CalculatePotentialGradient(Molecule[] moles)
    {
        float sig_6 = Mathf.Pow(sig, 6);
        Vector3 sum = new Vector3(0, 0, 0) ;
        for (int i = 0; i < moles.Length; i++)
        {
            if (moles[i] != this)
            {
                float sqrDistance = (Position - moles[i].Position).sqrMagnitude;
                float a = 0;
                try
                {
                    a = -24 * eps *
                    (2 * sig_6 * sig_6 / Mathf.Pow(sqrDistance, 7)
                     - sig_6 / Mathf.Pow(sqrDistance, 4));
                }
                catch
                {

                }

                Vector3 s = sum;
                sum += a * (moles[i].Position - Position);
                if(sum.x ==float.NaN)
                {
                    sum = s;
                }
            }
        }
        Debug.Log(sum);
        return sum;
    }

    float CalculatePursion(Molecule other)
    {
        float sqrDistance = Vector3.SqrMagnitude(Position - other.Position);

        float sig_6pow = (float)Mathf.Pow(sig, 6);
        float distance_6pow = (float)Mathf.Pow(sqrDistance, 3);
        return 4 * eps * (sig_6pow * sig_6pow / (distance_6pow * distance_6pow)
            - sig_6pow / distance_6pow);
    }

    public void Update(float deltaTime, Molecule[] moles, Vector3 boundaryCenter, Vector3 boundarySizeHalf)
    {
        Velocity = Velocity + CalculatePotentialGradient(moles).normalized;
        Position = Position + Velocity * deltaTime;

        float x = Position.x;
        if (x < boundaryCenter.x - boundarySizeHalf.x)
        {
            x = boundaryCenter.x - boundarySizeHalf.x;
            Velocity.x = -Velocity.x;
        }
        else if (x > boundaryCenter.x + boundarySizeHalf.x)
        {
            x = boundaryCenter.x + boundarySizeHalf.x;
            Velocity.x = -Velocity.x;
        }

        float y = Position.y;
        if (y < boundaryCenter.y - boundarySizeHalf.y)
        {
            y = boundaryCenter.y - boundarySizeHalf.y;
            Velocity.y = -Velocity.y;
        }
        else if (y > boundaryCenter.y + boundarySizeHalf.y)
        {
            y = boundaryCenter.y + boundarySizeHalf.y;
            Velocity.y = -Velocity.y;
        }

        float z = Position.z;
        if (z < boundaryCenter.z - boundarySizeHalf.z)
        {
            z = boundaryCenter.z - boundarySizeHalf.z;
            Velocity.z = -Velocity.z;
        }
        else if (z > boundaryCenter.z + boundarySizeHalf.z)
        { 
            z = boundaryCenter.z + boundarySizeHalf.z;
            Velocity.z = -Velocity.z;
        }

        Position = new Vector3(x, y, z);
    }
}
