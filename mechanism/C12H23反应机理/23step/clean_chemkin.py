# -*- coding: utf-8 -*-
import re

def strip_invalid_ford(inp_path, out_path_remove, out_path_threebody):
    with open(inp_path, 'r') as f:
        lines = f.readlines()

    out_lines_remove = []      # Version 1: Remove invalid FORD
    out_lines_threebody = []   # Version 2: Convert to three-body reactions
    net_stoich = {}           # 当前反应的净化学计量字典
    in_reactions = False
    current_reaction_line = None

    # 匹配反应方程式行 (包含动力学参数) - 支持 => 和 <=>
    reac_re = re.compile(r'^\s*([^=<\s]+)\s*<?=>\s*([^\s]+)\s+(.+)$')
    # 匹配 FORD 行并捕获物种名
    ford_re = re.compile(r'^\s*FORD/(\w+)\s')

    for L in lines:
        m = reac_re.match(L)
        if m:
            # 进入一个新反应，重置 net_stoich
            in_reactions = True
            net_stoich.clear()
            current_reaction_line = L

            lhs, rhs, kinetics = m.group(1), m.group(2), m.group(3)
            # 解析反应物和生成物
            def parse_side(side_str, sign):
                # sign: +1 for LHS, -1 for RHS
                parts = side_str.split('+')
                for term in parts:
                    term = term.strip()
                    # 提取前面的系数（可能为空，默认为1）
                    m2 = re.match(r'^(\d*)([A-Za-z0-9_]+)$', term)
                    if not m2:
                        continue
                    coef = int(m2.group(1)) if m2.group(1) else 1
                    sp   = m2.group(2)
                    net_stoich[sp] = net_stoich.get(sp, 0) + sign * coef

            parse_side(lhs, +1)
            parse_side(rhs, -1)

            out_lines_remove.append(L)
            out_lines_threebody.append(L)
            continue

        # 如果出到 REACTIONS 段外，直接保留
        if not in_reactions:
            out_lines_remove.append(L)
            out_lines_threebody.append(L)
            continue

        # 如果遇到下一段（如空行、END、元素或 SPECIES 段），退出 reactions 段
        if L.strip() == '' or re.match(r'^(END|ELEMENTS|SPECIES)', L):
            in_reactions = False
            out_lines_remove.append(L)
            out_lines_threebody.append(L)
            continue

        # 如果是 FORD 行
        m2 = ford_re.match(L)
        if m2:
            sp = m2.group(1)
            # Version 1: 只有当 net_stoich[sp] > 0（即真正被消耗）才保留
            if net_stoich.get(sp, 0) > 0:
                out_lines_remove.append(L)
                out_lines_threebody.append(L)
            else:
                # Version 2: Convert invalid FORD to three-body reaction
                # For threebody version, modify the reaction equation to include (+M)
                if current_reaction_line:
                    # Find the last reaction line in threebody output and modify it
                    for i in range(len(out_lines_threebody) - 1, -1, -1):
                        if reac_re.match(out_lines_threebody[i]):
                            # Add (+M) to the reaction if not already present
                            reaction_line = out_lines_threebody[i].strip()
                            if '(+M)' not in reaction_line:
                                # Split reaction equation and kinetics
                                m_temp = reac_re.match(reaction_line)
                                if m_temp:
                                    lhs_temp, rhs_temp, kinetics_temp = m_temp.group(1), m_temp.group(2), m_temp.group(3)
                                    # Add (+M) to both sides of reaction equation only
                                    new_reaction = lhs_temp.strip() + '(+M) => ' + rhs_temp.strip() + '(+M)   ' + kinetics_temp
                                    out_lines_threebody[i] = new_reaction + '\n'
                            break
                # Add the third-body efficiency line instead of FORD
                efficiency_line = "    {}/1.0/\n".format(sp)
                out_lines_threebody.append(efficiency_line)
            continue

        # 其他行一律保留
        out_lines_remove.append(L)
        out_lines_threebody.append(L)

    # 写出清理后的文件
    with open(out_path_remove, 'w') as f:
        f.writelines(out_lines_remove)
    
    with open(out_path_threebody, 'w') as f:
        f.writelines(out_lines_threebody)

if __name__ == '__main__':
    inp  = 'jeta-23steps.inp'
    outp_remove = 'jeta-23steps_clean_remove.inp'
    outp_threebody = 'jeta-23steps_clean_threebody.inp'
    strip_invalid_ford(inp, outp_remove, outp_threebody)
    print("Version 1 (remove invalid FORD): {}".format(outp_remove))
    print("Version 2 (convert to three-body): {}".format(outp_threebody))