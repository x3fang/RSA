#include "RSA.h"
ulli hcf_(ulli a, ulli b)
{
    return b == 0 ? a : hcf_(b, a % b);
}
inline ulli hcf(ulli a, ulli b)
{
    return a >= b ? hcf_(a, b) : hcf_(b, a);
}

inline ulli lcm(ulli a, ulli b)
{
    return (a * b) / hcf(a, b);
}
bool millerRabin(ulli n, ulli iterations); // Miller-Rabin 质数测试函数
bool isPrime(ulli n)                       // O(sqrt(n))
{
    if (n > MAX_RANGE)
    {
        return millerRabin(n, millerRabinTimes);
    }
    if (n <= 3)
        return n > 1;
    if (n % 6 != 1 && n % 6 != 5)
        return false;

    auto Sqrt = sqrt(n);
    for (ulli i = 5; i <= Sqrt; i += 6)
        if (n % i == 0 || n % (i + 2) == 0)
            return false;
    return true;
}
ulli quik_power(ulli a, ulli b, ulli n)
{
    ulli result = 1;
    while (b > 0)
    {
        if (b & 1)
            result = result * a % n;
        a = a * a % n;
        b >>= 1;
    }
    if (n == 0)
        throw n;
    return result % n;
}
void CalPThread(ulli range1, ulli range2, ulli &p_T);
ulli Rsa::CalNextP()
{
    if (range1 <= 0 || range2 <= 0)
        return 0;
    if (abs(max(range1, max(p, q) + 1) - range2) <= MAX_RANGE)
        for (ulli i = max(range1, max(p, q) + 1); i <= range2; i++)
        {
            if (i % 6 != 1 && i % 6 != 5)
                continue;
            if (isPrime(i))
            {
                p = i;
                return i;
            }
        }
    else
    {
        std::vector<std::thread> t;
        ulli min_p = range1 * range2 + 2;
        ulli p_T[MAX_THREADS + 2];
        ulli range1t = max(range1, max(p, q) + 1);
        ulli addNum = (range2 - range1) / MAX_THREADS;
        ulli range2t = range1 + addNum;
        for (int i = 1; i < MAX_THREADS; i++)
        {
            t.push_back(std::thread(CalPThread,
                                    range1t,
                                    range2t,
                                    std::ref(p_T[i])));
            range1t = range2t + 1;
            range2t = range1t + addNum;
        }
        for (auto &th : t)
            th.join();
        for (int i = 1; i < MAX_THREADS; i++)
            if (p_T[i] > 1 && p_T[i] < min_p)
                min_p = p_T[i];
        p = (min_p != (range1 * range2 + 2) ? min_p : 0);
        return p;
    }
    return 0;
}
void CalPThread(ulli range1, ulli range2, ulli &p_T)
{
    for (ulli i = range1; i <= range2; i++)
    {
        if (i % 6 != 1 && i % 6 != 5)
            continue;
        if (isPrime(i))
        {
            p_T = i;
            return;
        }
    }
    p_T = 0;
}
ulli Rsa::CalNextQ()
{
    if (range1 <= 0 || range2 <= 0)
        return 0;
    if (abs(max(range1, max(p, q) + 1) - range2) <= MAX_RANGE)
        for (ulli i = max(range1, max(p, q) + 1); i <= range2; i++)
        {
            if (i % 6 != 1 && i % 6 != 5)
                continue;
            if (isPrime(i) && i != p)
            {
                q = i;
                return i;
            }
        }
    else
    {
        std::vector<std::thread> t;
        ulli min_q = range1 * range2 + 2;
        ulli q_T[MAX_THREADS + 2];
        ulli range1t = max(range1, max(p, q) + 1);
        ulli addNum = (range2 - range1) / MAX_THREADS;
        ulli range2t = range1 + addNum;
        for (int i = 1; i < MAX_THREADS; i++)
        {
            t.push_back(std::thread(CalPThread,
                                    range1t,
                                    range2t,
                                    std::ref(q_T[i])));
            range1t = range2t + 1;
            range2t = range1t + addNum;
        }
        for (auto &th : t)
            th.join();
        for (int i = 1; i < MAX_THREADS; i++)
            if (q_T[i] > 1 && q_T[i] < min_q && q_T[i] != p)
                min_q = q_T[i];
        q = (min_q != (range1 * range2 + 2) ? min_q : 0);
        return q;
    }
    return 0;
}
void CalEThread(ulli &phi, ulli range1, ulli range2, ulli &e_T)
{
    for (ulli i = range1; i <= range2; i++)
    {
        if (hcf(i, phi) == 1)
        {
            e_T = i;
            return;
        }
    }
    e_T = 0;
}
ulli Rsa::CalNextE()
{
    if (p <= 0 || q <= 0)
        return 0;
    phi = (p - 1) * (q - 1);
    n = p * q;
    if (phi <= MAX_RANGE)
        for (ulli i = 2; i < phi; i++)
        {
            if (hcf(i, phi) == 1)
            {
                e = i;
                return i;
            }
        }
    else
    {
        std::vector<std::thread> t;
        ulli min_e = range1 * range2 + 2;
        ulli e_T[MAX_THREADS + 2];
        ulli range1t = max(range1, max(p, q) + 1);
        ulli addNum = (range2 - range1) / MAX_THREADS;
        ulli range2t = range1t + addNum;
        for (int i = 1; i < MAX_THREADS; i++)
        {
            t.push_back(std::thread(CalEThread,
                                    std::ref(phi),
                                    range1t,
                                    range2t + (i == MAX_THREADS - 1 ? -1 : 0),
                                    std::ref(e_T[i])));
            range1t = range2t + 1;
            range2t = range1t + addNum;
        }
        for (auto &th : t)
            th.join();
        for (int i = 1; i < MAX_THREADS; i++)
            if (e_T[i] > 1 && e_T[i] < min_e)
                min_e = e_T[i];
        e = (min_e != (range1 * range2 + 2) ? min_e : 0);
        return e;
    }
    return 0;
}
ulli calc_d(ulli e, ulli n)
{
    if (e == 1 || n == 1)
        return 1;
    ulli d1 = calc_d(e % n, n % e);
    return d1 - n / e * ((1 - e * d1) / (n % e));
}
ulli Rsa::CalNextD()
{
    if (e <= 0 || phi <= 0)
        return 0;
    d = calc_d(e, phi);
    return d;
}

bool Rsa::CalPQED()
{
    if (CalNextP())
        if (CalNextQ())
            if (CalNextE())
                if (CalNextD())
                    return true;
    return false;
}
std::vector<ulli> Encrypt(ulli e, ulli n, ulli plainText)
{
    std::vector<ulli> cipherText;
    if (plainText > n)
    {
        ulli temp = -1;

        std::string str = to_string(plainText);
        for (auto i : str)
        {
            if ((temp * 10 + (i - '0')) <= n && temp > 0)
            {
                temp = temp * 10 + (i - '0');
            }
            else if (temp == -1)
            {
                temp = (i - '0');
            }
            else
            {
                cipherText.push_back(temp);
                temp = (i - '0');
            }
        }
        cipherText.push_back(temp);
        for (auto &i : cipherText)
            i = quik_power(i, e, n);
        std::cout << cipherText.size() << std::endl;
    }
    return cipherText;
}
ulli Decrypt(ulli d, ulli n, std::vector<ulli> cipherText)
{
    std::string str;
    for (auto i : cipherText)
        str += to_string(quik_power(i, d, n));
    return ulli(str);
}
// Miller-Rabin 质数测试函数，参数为 n（待检测数）和 iterations（测试的迭代次数）
bool millerRabin(ulli n, ulli iterations)
{
    if (n < 2) // 如果 n 小于 2，则不是素数
        return false;

    // 将 n-1 表示为 2^s * d 的形式，r 是 n-1，不断除以2
    ulli r = n - 1;
    while (r % 2 == 0) // r 是偶数时，继续除以 2
        r /= 2;

    // 进行指定次数的迭代测试
    for (ulli i = 0; i < iterations; i++)
    {
        // 随机选择 a，范围在 [2, n-2] 之间
        ulli a = 2 + rand() % (n - 4);

        // 计算 x = a^r % n
        ulli x = quik_power(a, r, n);
        ulli temp = r;

        // 如果 x == 1 或 x == n-1，则继续下一次迭代
        if (x == 1 || x == n - 1)
            continue;

        // 检查是否为合数
        bool composite = true;
        // 重复平方，检查 x 是否等于 n-1
        while (temp != n - 1)
        {
            x = (x * x) % n; // 更新 x = x^2 % n
            temp *= 2;       // temp = temp * 2

            // 如果 x == 1，则 n 不是素数
            if (x == 1)
                return false;

            // 如果 x == n-1，则 n 可能是素数
            if (x == n - 1)
            {
                composite = false; // 不是合数，继续进行测试
                break;
            }
        }

        // 如果所有测试都未找到 n 是素数的证据，则 n 是合数
        if (composite)
            return false;
    }

    // 如果所有测试通过，返回 n 是素数
    return true;
}
