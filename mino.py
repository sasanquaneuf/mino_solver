from datetime import datetime
from typing import List, Tuple, Set, Dict
from pprint import pprint
from json import loads

debug = False


def normalize(mino: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    min_x = min([x for x, _ in mino])
    min_y = min([y for _, y in mino])
    return {(x - min_x, y - min_y) for x, y in mino}


def rotate(mino: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    return {(-y, x) for x, y in mino}


def reverse(mino: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    return {(-x, y) for x, y in mino}


def semi_standardize(mino: Set[Tuple[int, int]]) -> List[Tuple[Set[Tuple[int, int]], Tuple[int, int]]]:
    semi_standard_mino_list: List[Tuple[Set[Tuple[int, int]], Tuple[int, int]]] = []
    for wx, wy in mino:
        if (wx - 1, wy) not in mino and (wx, wy - 1) not in mino:
            semi_standard_mino_list.append((
                {(x - wx, y - wy) for x, y in mino},
                (wx, wy)))
    return semi_standard_mino_list


def get_key(mino: Set[Tuple[int, int]]) -> str:
    return str(sorted(list(mino), key=lambda a: a[1] * 10000 + a[0]))


def str_mino(mino: Set[Tuple[int, int]]) -> str:
    # normalized
    output_str = ''
    min_x = prev_x = min([x for x, _ in mino]) - 1
    prev_y = min([y for _, y in mino])
    for x, y in sorted(list(mino), key=lambda a: a[1] * 10000 + a[0]):
        if prev_y != y:
            output_str += '\n'
            prev_x = min_x
        output_str += '　' * (x - prev_x - 1) + '箱'
        prev_x = x
        prev_y = y
    return output_str


def append_mino(
        working_mino: Set[Tuple[int, int]], count: int,
        mino_representative_list: List[str], key_to_representative: Dict[str, str],
        key_to_mino: Dict[str, Set[Tuple[int, int]]]) -> None:
    if count == 0:
        mino = normalize(working_mino)
        representative_key = get_key(mino)
        if representative_key not in key_to_representative:
            if debug:
                print('---- new mino ----')
                print(mino)
                print(str_mino(mino))
                print('------------------')
            mino_representative_list.append(representative_key)
            for _ in range(2):
                mino = reverse(mino)
                for __ in range(4):
                    mino = normalize(rotate(mino))
                    mino_key = get_key(mino)
                    if mino_key not in key_to_representative:
                        key_to_representative[mino_key] = representative_key
                        key_to_mino[mino_key] = mino
        return
    for wx, wy in working_mino:
        for dx, dy in [(1, 0), (0, 1), (0, -1)]:
            if (wx + dx, wy + dy) not in working_mino:
                mino = working_mino.copy()
                mino.add((wx + dx, wy + dy))
                append_mino(
                    mino, count - 1, mino_representative_list,
                    key_to_representative, key_to_mino)


def append_ss_mino(
        key_to_representative: Dict[str, str],
        key_to_mino: Dict[str, Set[Tuple[int, int]]],
        ss_key_to_representative: Dict[str, str],
        ss_key_to_normal: Dict[str, str],
        ss_key_to_mino: Dict[str, Set[Tuple[int, int]]],
        ss_key_to_diff: Dict[str, Tuple[int, int]],
        representative_to_ss_mino_list: Dict[str, List[Tuple[Set[Tuple[int, int]], Tuple[int, int], str]]]
        ) -> None:
    for key, mino in key_to_mino.items():
        if key_to_representative[key] not in representative_to_ss_mino_list:
            representative_to_ss_mino_list[key_to_representative[key]] = []
        for semi_standard_mino, diff in semi_standardize(mino):
            ss_key = get_key(semi_standard_mino)
            ss_key_to_representative[ss_key] = key_to_representative[key]
            ss_key_to_normal[ss_key] = key
            ss_key_to_mino[ss_key] = semi_standard_mino
            ss_key_to_diff[ss_key] = diff
            representative_to_ss_mino_list[key_to_representative[key]].append((semi_standard_mino, diff, key))


def solve(
        board: List[List[int]],
        count: int, placed_normal_list: List[int],
        calculated_str_set: Set[str], solved_str_set: Set[str],
        left_top_filled_cells: List[Tuple[int, int]],
        representative_left_dict: Dict[str, int],
        maximum_mino_count: int,
        mino_representative_list: List[str],
        representative_to_ss_mino_list: Dict[str, List[Tuple[Set[Tuple[int, int]], Tuple[int, int], str]]],
        representative_to_num: Dict[str, int],
        key_to_num: Dict[str, int]) -> None:
    if str(placed_normal_list) in calculated_str_set:
        return
    calculated_str_set.add(str(placed_normal_list))
    if count == 0:
        exact_placed_dict: Dict[int, List[Tuple[int, int]]] = {}
        for i in range(len(key_to_num)):
            for j in range(maximum_mino_count):
                temp_position = (i * maximum_mino_count + j) * 2
                if placed_normal_list[temp_position] != -1:
                    if i not in exact_placed_dict:
                        exact_placed_dict[i] = []
                    exact_placed_dict[i].append((
                        placed_normal_list[temp_position],
                        placed_normal_list[temp_position + 1]
                    ))
        for value in exact_placed_dict.values():
            value.sort()
        exact_placed_str = str(sorted(list(exact_placed_dict.items())))
        if exact_placed_str not in solved_str_set:
            solved_str_set.add(exact_placed_str)
            if debug:
                print('---- solved ----')
                print(datetime.now())
                print(exact_placed_str)
                print(len(solved_str_set))
                pprint(board)
                print('----------------')
        return
    # 埋めたセルの右または下を見る
    for cell_x, cell_y in left_top_filled_cells:
        for representative in mino_representative_list:
            if representative_left_dict[representative] > 0:
                # はまるミノを探す
                representative_num = representative_to_num[representative]
                for ss_mino, diff, key in representative_to_ss_mino_list[representative]:
                    filled = False
                    for dx, dy in ss_mino:
                        if board[cell_y + dy][cell_x + dx] != 0:
                            filled = True
                    if not filled:
                        next_left_top_filled_cells = list(left_top_filled_cells)
                        for dx, dy in ss_mino:
                            board[cell_y + dy][cell_x + dx] = representative_num
                            if (cell_x + dx, cell_y + dy) in next_left_top_filled_cells:
                                next_left_top_filled_cells.remove((cell_x + dx, cell_y + dy))
                        for dx, dy in ss_mino:
                            # 埋めたセルの右上が埋まって右が空、または左下が埋まって下が空
                            temp_x = cell_x + dx
                            temp_y = cell_y + dy
                            if board[temp_y][temp_x + 1] == 0 and board[temp_y - 1][temp_x + 1] != 0:
                                next_left_top_filled_cells.append((temp_x + 1, temp_y))
                            if board[temp_y + 1][temp_x] == 0 and board[temp_y + 1][temp_x - 1] != 0:
                                next_left_top_filled_cells.append((temp_x, temp_y + 1))
                        key_num = key_to_num[key]
                        i = -1
                        while (i := i + 1) < maximum_mino_count:
                            temp_position = (key_num * maximum_mino_count + i) * 2
                            if placed_normal_list[temp_position] == -1:
                                placed_normal_list[temp_position] = cell_x - diff[0]
                                placed_normal_list[temp_position + 1] = cell_y - diff[1]
                                i = maximum_mino_count
                        representative_left_dict[representative] -= 1
                        count -= len(ss_mino)
                        solve(
                            board, count, placed_normal_list, calculated_str_set, solved_str_set,
                            next_left_top_filled_cells, representative_left_dict, maximum_mino_count,
                            mino_representative_list, representative_to_ss_mino_list,
                            representative_to_num, key_to_num)
                        count += len(ss_mino)
                        representative_left_dict[representative] += 1
                        i = maximum_mino_count
                        while (i := i - 1) > -1:
                            temp_position = (key_num * maximum_mino_count + i) * 2
                            if placed_normal_list[temp_position] != -1:
                                placed_normal_list[temp_position] = -1
                                placed_normal_list[temp_position + 1] = -1
                                i = 0
                        for dx, dy in ss_mino:
                            board[cell_y + dy][cell_x + dx] = 0


def rotate_board(board: List[List[int]]) -> List[List[int]]:
    return [[cell for cell in row[::-1]] for row in board[::-1]]


def reverse_board(board: List[List[int]]) -> List[List[int]]:
    return [[cell for cell in row[::-1]] for row in board]


def solve_tetromino_puzzle(
        output_placement_data: bool = True, reduce_symmetric_data: bool = True, file_name: str = 'result.txt') -> int:
    mino_representative_list: List[str] = []
    key_to_representative: Dict[str, str] = {}
    key_to_mino: Dict[str, Set[Tuple[int, int]]] = {}
    append_mino(
        {(0, 0)}, 3, mino_representative_list, key_to_representative, key_to_mino)

    representative_to_num: Dict[str, int] = {
        key: index + 1 for index, key in enumerate(mino_representative_list)}
    key_to_num: Dict[str, int] = {
        key: index for index, key in enumerate(list(key_to_representative.keys()))}
    num_to_key: Dict[int, str] = {v: k for k, v in key_to_num.items()}

    ss_key_to_representative: Dict[str, str] = {}
    ss_key_to_normal: Dict[str, str] = {}
    ss_key_to_mino: Dict[str, Set[Tuple[int, int]]] = {}
    ss_key_to_diff: Dict[str, Tuple[int, int]] = {}
    representative_to_ss_mino_list: Dict[str, List[Tuple[Set[Tuple[int, int]], Tuple[int, int], str]]] = {}
    append_ss_mino(
        key_to_representative, key_to_mino,
        ss_key_to_representative, ss_key_to_normal, ss_key_to_mino,
        ss_key_to_diff, representative_to_ss_mino_list)

    if debug:
        print(mino_representative_list)
        print(key_to_representative)
        print(key_to_mino)

        for key, mino in ss_key_to_mino.items():
            print('----' * 5)
            print(str_mino(mino))
            print('is same as')
            print(str_mino(ss_key_to_mino[ss_key_to_representative[key]]))

    width = 8
    height = 5
    padding = 3
    board: List[List[int]] = [[-1] * padding + [0 for _ in range(width)] + [-1] * padding for __ in range(height)]
    board = [[-1] * (2 * padding + width)] * padding + board + [[-1] * (2 * padding + width)] * padding
    count = sum([len([cell for cell in row if cell == 0]) for row in board])
    left_top_filled_cells = [(padding, padding)]
    representative_left_dict = {representative: 2 for representative in mino_representative_list}
    maximum_mino_count = max(list(representative_left_dict.values()))
    # coordinate x and y -> 2, each mino's count -> maximum_mino_count
    placed_normal_list = [-1] * len(key_to_representative) * 2 * maximum_mino_count
    calculated_str_set: Set[str] = set()
    solved_str_set: Set[str] = set()

    if output_placement_data:
        solve(
            board,
            count,
            placed_normal_list,
            calculated_str_set,
            solved_str_set,
            left_top_filled_cells,
            representative_left_dict,
            maximum_mino_count,
            mino_representative_list,
            representative_to_ss_mino_list,
            representative_to_num,
            key_to_num,
            )

        with open(file_name, mode='w') as f:
            for solved_str in solved_str_set:
                f.write(solved_str + '\n')

    if reduce_symmetric_data:
        with open(file_name) as f:
            symmetric_free_list = []
            added_set = set()
            for line in f.readlines():
                placed_normal_position_list = loads(line.replace('(', '[').replace(')', ']'))
                board = [
                    [-1] * ((2 * padding + width) * 2 + 1) for _ in range((2 * padding + height) * 2 + 1)]
                for normal_num, position_list in placed_normal_position_list:
                    representative_num = representative_to_num[key_to_representative[num_to_key[normal_num]]]
                    mino = key_to_mino[num_to_key[normal_num]]
                    for x, y in position_list:
                        for dx, dy in mino:
                            board[(y + dy) * 2 + 1][(x + dx) * 2 + 1] = representative_num
                            if (dx - 1, dy) in mino:
                                board[(y + dy) * 2 + 1][(x + dx) * 2] = representative_num
                            if (dx, dy - 1) in mino:
                                board[(y + dy) * 2][(x + dx) * 2 + 1] = representative_num
                board_key = str(board)
                if board_key not in added_set:
                    added_set.add(board_key)
                    added_set.add(str(rotate_board(board)))
                    added_set.add(str(reverse_board(board)))
                    added_set.add(str(rotate_board(reverse_board(board))))
                    symmetric_free_list.append(board_key)
        return len(symmetric_free_list)
    else:
        return len(solved_str_set)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Solve tetromino puzzle.')
    parser.add_argument(
        '--debug', dest='debug', action='store_const',
        const=True, default=False,
        help='run as debug mode(default: off)')
    args = parser.parse_args()
    debug = args.debug
    print('This process may use 16GB of memory...')
    if not debug:
        print('If you want to watch detail, please set argument "--debug".')
    pattern_number = solve_tetromino_puzzle()
    print(pattern_number)
