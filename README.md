# Nonexistence Result of Lattice Tiling by Lee Balls
Computation Results of Article: On Lattice Tilings by Lee Balls: Q-Polynomials, Group Ring Equations, and a Decision Algorithm

### Maple代码示例

```maple
# 计算Lee球的函数
compute_lee_ball := proc(n, r)
    local ii, result;
    result := expand(add(2^ii * binomial(n, ii) * binomial(r, ii), ii = 0 .. r));
    return result;
end proc;
```

如果 `maple` 标签不被识别，尝试替代方案：
```javascript
// 使用JavaScript高亮作为替代
compute_lee_ball := proc(n, r)
    local ii, result;
    result := expand(add(2^ii * binomial(n, ii) * binomial(r, ii), ii = 0 .. r));
    return result;
end proc;
```

```text
# 纯文本展示，无高亮
compute_lee_ball := proc(n, r)
    local ii, result;
    ...
end proc;
```
